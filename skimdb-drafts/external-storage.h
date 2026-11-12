#pragma once

#include "refseq/dna.h"
#include "definitions.h"
#include "utilities.h"
#include "zipsort/zipsort.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <future>
#include <limits>
#include <memory>
#include <print>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#include <theta_sketch.hpp>
#include <theta_intersection.hpp>

#include <oneapi/tbb/info.h>
#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/concurrent_queue.h>

using select_type = kset_t::select_1_type;
using sketch_type = datasketches::update_theta_sketch;

/**
 * @brief Unified storage for intermediate kmer sets. Encapsulates POSIX API.
 */
class posix_kmer_set_file final {
  static constexpr std::size_t max_io_retries = 3;
  static constexpr std::string_view file_name = "kmers.tmp";

public:
  posix_kmer_set_file(const fs::path &dir)
      : path(dir / file_name),
        fd_(::open(path.c_str(), O_CREAT | O_RDWR | O_CLOEXEC | O_EXCL, S_IRUSR | S_IWUSR))
  {
    if (fd_ == -1) {
      throw std::system_error{errno, std::generic_category(), "open"};
    }
  }

  ~posix_kmer_set_file() {
    if (fd_ != -1) {
      ::close(fd_);
    }

    std::error_code ec;
    fs::remove(path, ec);

    if (ec) {
      std::println("Could not remove tmp file: {}", ec.message());
    }
  }

  posix_kmer_set_file(const posix_kmer_set_file&) = delete;
  auto operator=(const posix_kmer_set_file&) -> posix_kmer_set_file& = delete;
  posix_kmer_set_file(posix_kmer_set_file&&) = delete;
  auto operator=(posix_kmer_set_file&&) -> posix_kmer_set_file& = delete;

  /**
   * @brief Synchronously write an indexed buffer of kmer sets to external storage.
   * Throws on write failures.
   */
  void sync_write(indexed_stream& is) {
    const std::size_t base = eof_;
    std::size_t retries = 0;
    std::string_view raw = is.str.view();
    
    while (!raw.empty()) {
      const auto n = ::pwrite(fd_, raw.data(), raw.size(), static_cast<off_t>(eof_));
      
      if (n > 0) [[likely]] {
        raw.remove_prefix(n);
        eof_ += static_cast<std::size_t>(n);
      }
      else if (n == 0 && ++retries > max_io_retries) {
        throw std::runtime_error{std::format(
            "fatal: zero write, but expected {} bytes.", raw.size())};
      }
      else if (n < 0 && errno != EINTR && ++retries > max_io_retries) {
        throw std::system_error{errno, std::generic_category(), "pwrite"};
      }
    }

    for (auto meta : is.idx) {
      meta.get().shift_offset(base);
    }

    is.str.str("");
    is.str.clear();
    is.idx.clear();
  }

  /**
   * @brief Synchronously load kmer set from external storage.
   * Throws on read failures.
   */
  void sync_read(char* dst, std::size_t len, std::uint64_t off) const {
    std::size_t fetched = 0;
    std::size_t retries = 0;

    while (fetched < len) {
      const auto rc = ::pread(
          fd_, dst + fetched, len - fetched, static_cast<off_t>(off + fetched)
      );

      if (rc > 0) [[likely]] {
        fetched += static_cast<std::size_t>(rc);
      }
      else if (rc < 0 && errno == EINTR) {
        continue;
      }
      else if (rc == 0) {
        throw std::runtime_error{std::format(
            "fatal: expected {} bytes, but got {}.", len, fetched)};
      }
      else if (++retries > max_io_retries) {
        throw std::system_error{errno, std::generic_category(), "pread"};
      }
    }
  }

private:
  fs::path path;
  std::uint64_t eof_ = 0;
  int fd_ = -1;
};

/**
 * @brief Operates a fixed-size window of current kmer sets in memory.
 * Automates sketching and replacement for TSP sorting.
 */
class sequence_ring final {
private:
  /**
   * @brief Computes relative ring buffer position from an absolute index.
   */
  auto slot_for(std::uint32_t idx) const -> std::uint32_t {
    if (idx <= pointer_ || idx >= view_.size()) [[unlikely]] {
      throw std::runtime_error{std::format(
          "fatal: index ({}) outside range ({} - {})",
          idx, pointer_, view_.size() - 1)};
    }

    const std::size_t rel = idx - pointer_;

    if (rel >= size_) [[unlikely]] {
      throw std::runtime_error{std::format(
          "fatal: index ({}) exceeds ring size ({})", idx, size_)};
    }

    return static_cast<std::uint32_t>(
        (static_cast<std::size_t>(start_) + rel) % size_);
  }

  /**
   * @brief Executes a callback after IO completion. Can throw a stored
   * IO task exception or an error if IO states do not match.
   */
  template <typename F>
  auto after_io(std::uint32_t idx, F&& def) const -> std::uint32_t {
    // Contract on future completion is strict; spurious wakeups
    // have already been accounted for by the implementation.
    drdy_.get();

    const kview::io_status st = view_[idx].state();

    if (st != kview::Present) [[unlikely]] {
      throw std::runtime_error{std::format(
          "fatal: unparked but kmer set is not marked ready ({})",
          idx, kview::io_to_str(st))};
    }

    return def();
  }

  /**
   * @brief Synchronously reconstructs a kmer set from external storage
   * and populates its sketch.
   */
  void prepare_kmer_set(std::uint32_t idx) {
    const std::size_t len = view_[idx].size;
    const std::size_t bits_set = view_[idx].card;
    const std::uint32_t slot = slot_for(idx);

    sketches_[slot].reset();

    auto& set = sets_[iota_[slot]];

    if (len == 0) [[unlikely]] {
      sdsl::util::clear(set);
      return;
    }

    std::string raw(len, '\0');
    file_.sync_read(raw.data(), len, view_[idx].offset());

    std::istringstream in{std::move(raw), std::ios::in | std::ios::binary};
    {
      cereal::BinaryInputArchive archive{in};
      archive(set);
    }

    const select_type select{&set};

    for (std::size_t i = 1; i <= bits_set; ++i) {
      sketches_[slot].update(select(i));
    }
  }

public:
  sequence_ring() = delete;

  sequence_ring(std::size_t s,
                std::vector<kview>& v,
                const posix_kmer_set_file& f)
  // No allocation or setup -- those must only be done when necessary.
      : size_(s + 1), view_(v), file_(f) {}

  ~sequence_ring() {
    if (!drdy_.valid()) {
      return;
    }

    try {
      drdy_.get();
    } catch (const std::exception& e) {
      // We are finishing anyway, and one does not throw from a destructor :).
      std::println("Caught an exception from async routine: {}", e.what());
    }
  }

  sequence_ring(const sequence_ring&) = delete;
  auto operator=(const sequence_ring&) -> sequence_ring& = delete;
  sequence_ring(sequence_ring&&) = delete;
  auto operator=(sequence_ring&&) -> sequence_ring& = delete;

  /**
   * @brief Allows to suspend processing inside the ring buffer until
   * signaled using the provided handle.
   */
  auto defer() const noexcept -> std::promise<void> {
    std::promise<void> handle;
    drdy_ = handle.get_future().share();
    return handle;
  }

  /**
   * @brief Allocates ring buffers, adjust metadata, and initiates
   * an async transfer to fill the entire buffer.
   */
  void initialize() {
    sets_ = std::make_unique<kset_t[]>(size_);
    iota_ = std::make_unique<std::uint32_t[]>(size_);
    sketches_.reserve(size_);

    sketch_type::builder bd;

    // 2^13 = 8-kb register storage per window item.
    bd.set_lg_k(13);

    for (std::size_t i = 0; i < size_; ++i) {
      sketches_.emplace_back(bd.build());
      // Set set idx to avoid swapping sets.
      iota_[i] = i;
    }

    const std::uint32_t count = std::min(
      size_, static_cast<std::uint32_t>(view_.size())
    );

    if (drdy_.valid()) {
      drdy_.get();
    }

    for (std::size_t idx = 0; idx < count; ++idx) {
      view_[idx].set_state(kview::Pending);
    }

    drdy_ = std::async(
        std::launch::async,
        [this, count] {
          for (std::size_t idx = 0; idx < count; ++idx) {
            prepare_kmer_set(idx);
            view_[idx].set_state(kview::Present);
          }
        }
    ).share();
  }

  /**
   * @brief Computes cardinality of intersection between base set and a
   * set at provided index. Executes immediately if a set is present
   * or suspends current thread (!) until IO completes.
   */
  auto intersect(std::uint32_t idx) const -> std::uint32_t {
    std::uint32_t rb_idx = slot_for(idx);

    datasketches::theta_intersection op;
    op.update(sketches_[start_]);

    const auto is = [rb_idx, &op, this]() -> std::uint32_t {
      op.update(sketches_[rb_idx]);
      return static_cast<std::uint32_t>(op.get_result().get_estimate());
    };

    const kview::io_status st = view_[idx].state();

    switch (st) {
      case kview::Present:
        return is();
      case kview::Pending:
        return after_io(idx, std::move(is));
      case kview::External:
      default:
        throw std::runtime_error{std::format(
            "fatal: unexpected IO state: {}", kview::io_to_str(st))};
    }
  }

  /**
   * @brief Advances ring buffer and swaps next base set with the best candidate.
   * Initiates IO for the next set. Returns evicted set and its cardinality.
   */
  auto advance(std::uint32_t best_pos) -> std::pair<kset_t, std::size_t> {
    const std::uint32_t rb_best = slot_for(best_pos);

    kset_t evicted = std::move(sets_[iota_[start_]]);

    const std::size_t card = view_[best_pos].card;
    const std::size_t next = pointer_ + size_;

    start_ = (start_ + 1) % size_;
    ++pointer_;

    std::swap(iota_[start_], iota_[rb_best]);
    std::swap(sketches_[start_], sketches_[rb_best]);

    if (next < view_.size()) {
      view_[next].set_state(kview::Pending);

      drdy_ = std::async(
          std::launch::async,
          [this, next] {
            prepare_kmer_set(next);
            view_[next].set_state(kview::Present);
          }
      ).share();
    }
    else {
      drdy_ = {};
    }

    return {std::move(evicted), card};
  }

private:
  std::uint32_t start_ = 0;
  const std::uint32_t size_;
  std::uint64_t pointer_ = 0;

  std::unique_ptr<kset_t[]> sets_{};
  std::unique_ptr<std::uint32_t[]> iota_{};
  std::vector<sketch_type> sketches_{};

  std::vector<kview>& view_;
  const posix_kmer_set_file& file_;

  mutable std::shared_future<void> drdy_;
};

/**
 * @brief An abstraction over intermediate per-sequence kmer storage.
 *
 * Combines kmer file, sequence ring, writer event loop, and writer routine.
 * Catches and propagates multi-threaded exceptions within main stack frame.
 */
class external_storage final {
private:
  /**
   * @brief Event loop used for populating kmer file. Upon exiting,
   * signals sequence ring that it can begin computing.
   */
  void io_writer_event_loop(std::promise<void> handle) {
    int ticket = -1;

    for (;;) {
      // Worker blocks on empty queue.
      writeq_.pop(ticket);
      if (ticket == -1) {
        break;
      }

      file_.sync_write(writebuf_[ticket]);
      tickets_.give(ticket);
      ticket = -1;
    }

    // Flush all stale caches.
    for (int i = 0; i < vacant_; ++i) {
      const auto& set = writebuf_[i].idx;

      if (!set.empty()) {
        file_.sync_write(writebuf_[i]);
        tickets_.give(i);
      }
    }

    // Wake the ring to let it initialize.
    handle.set_value();
    writebuf_ = nullptr;
  }

public:
  external_storage() = delete;

  explicit external_storage(std::vector<kview>& kmers,
                            const fs::path &dir,
                            std::size_t win_size,
                            std::size_t ram_limit)            
      : file_(dir),
        ring_(win_size, kmers, file_),
        // We split available memory into blocks whose count is twice the expected
        // amount of workers in TBB-backed parallel algorithms, with a hope that
        // when a producer thread submits its local stream to the io worker, it can
        // grab a second block and start populating it immediately. Queue capacity
        // is, therefore, set to maximum number of blocks that can be submitted.
        vacant_(tbb::info::default_concurrency() * 2),
        worker_quota_(ram_limit / sizeof(kmer_t)),
        tickets_(vacant_),
        writebuf_(std::make_unique<indexed_stream[]>(vacant_)),
        george_rr_martin([this] { io_writer_event_loop(ring_.defer()); })
  {
    writeq_.set_capacity(vacant_);
  }

  ~external_storage() {
    george_rr_martin.join();
  }

  external_storage(const external_storage&) = delete;
  auto operator=(const external_storage&) -> external_storage& = delete;
  external_storage(external_storage&&) = delete;
  auto operator=(external_storage&&) -> external_storage& = delete;

  /**
   * @brief Tells George R.R. Martin to wake up and quit.
   */
  void signal_write_completion() {
    writeq_.push(-1);
  }

  /**
   * @brief Starts async transfer of initial candidates, in amount
   * equal to window size.
   */
  void fetch_starting_window() {
    ring_.initialize();
  }

  /**
    * @brief A wrapper around ring intersection to save an exception until
    * it is checked for within the main stack frame. Because exceptions are
    * only thrown for programming errors, main frame should safely terminate.
    */
  auto base_intersect(std::uint32_t j) {
    std::uint32_t res = 0;

    try {
      res = ring_.intersect(j);
    } catch (...) {
      cancellation = std::current_exception();
    }

    return res;
  }

  /**
   * @brief Consumes the kmer set and updates all database rows
   * for which that set contains the kmer.
   */
  void append_column(fakedb& db, std::uint32_t seq) {
    if (cancellation != nullptr) {
      // Rethrow inside main stack frame.
      std::rethrow_exception(cancellation);
    }

    const auto [kset, card] = ring_.advance(seq);

    const select_type select{&kset};
    
    try {
      tbb::parallel_for(
        tbb::blocked_range<kset_t::size_type>(1, card + 1),
        [&](const auto& range) {
          for (auto j = range.begin(); j != range.end(); ++j) {
            const kmer_t kmer = select(j);
            db.data[db.hash.find(kmer).value()].push_back(seq);
          }
        }
      );
    } catch (...) {
      // Rethrow within main stack frame
      std::rethrow_exception(std::current_exception());
    }
  }

  /**
    * @brief Serializes kmer set, updates its metadata, and places
    * the buffer into a local queue. Can park the thread if low on RAM.
    */
  void write(kset_t ks, std::size_t card, kview& meta) {
    thread_local int key = -1;

    if (key == -1) {
      // Can block when no vacant blocks are available.
      // Here we are doing this to avoid workers going OOM.
      key = tickets_.take();
    }

    const auto start = writebuf_[key].str.tellp();

    {
      cereal::BinaryOutputArchive archive{writebuf_[key].str};
      archive(ks);
      sdsl::util::clear(ks);
    }

    const auto end = writebuf_[key].str.tellp();

    if (!writebuf_[key].str || start == -1 || end <= start) [[unlikely]] {
      throw std::runtime_error{"fatal: write serialization failed"};
    }

    const std::uint64_t len = end - start;

    constexpr std::size_t max_u32 =
        static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max());
        
    if (len > std::min(worker_quota_, max_u32)) [[unlikely]] {
      throw std::runtime_error{std::format(
          "fatal: kmer set size ({}) exceeds quota ({})",
          len, worker_quota_)};
    }

    meta = kview{static_cast<std::uint64_t>(start), len, card};

    writebuf_[key].idx.push_back(std::ref(meta));
    
    if (end >= worker_quota_) {
      // Should not park because if there is a vacant block
      // to submit, then queue is not at capacity.
      writeq_.push(key);
      key = -1;
    }
  }

private:
  posix_kmer_set_file file_;
  sequence_ring ring_;

  const std::uint32_t vacant_;
  std::size_t worker_quota_ = 0;
  self_checkout tickets_;

  std::unique_ptr<indexed_stream[]> writebuf_ = nullptr;
  tbb::concurrent_bounded_queue<int> writeq_;
  std::thread george_rr_martin;

  std::exception_ptr cancellation = nullptr;
};


static inline void build_index(
    fakedb& db,
    external_storage& kmers,
    std::vector<kview>& kv,
    std::future<void> master_set_builder,
    std::uint32_t win_size = 0
) {
  std::uint32_t n = kv.size();

  if (n < 3 || kv.size() > std::numeric_limits<std::uint32_t>::max()) {
    throw std::runtime_error{std::format("fatal: too many kmers: {}", kv.size())};
  }

  win_size = win_size == 0 ? n - 1 : std::min(win_size, n - 1);

  zipsort(std::execution::par, kv, db.sequences, [](const auto& a, const auto& b) {
    return a.first.card > b.first.card;
  });

  kmers.fetch_starting_window();

  std::vector<std::uint32_t> idx(win_size);

  master_set_builder.get();

  for (std::uint32_t i = 0, end = n - 2; i < end; ++i) {
    std::uint32_t const ref_idx = i + 1;

    idx.resize(std::min(win_size, n - ref_idx));
    std::iota(idx.begin(), idx.end(), ref_idx);

    auto [best_card, best_pos] = std::transform_reduce(
        std::execution::par,
        idx.begin(),
        idx.end(),
        std::pair(std::uint32_t{0}, ref_idx),
        [](const auto a, const auto b) {
          return a.first >= b.first ? a : b;
        },
        [&](std::uint32_t j) -> std::pair<std::uint32_t, std::uint32_t> {
          return {kmers.base_intersect(j), j};
        });

    std::swap(kv[ref_idx], kv[best_pos]);
    std::swap(db.sequences[ref_idx], db.sequences[best_pos]);

    kmers.append_column(db, i);
  }

  // Real database would need to be compressed here.
}