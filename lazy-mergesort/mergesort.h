#pragma once

#include "refseq/dna.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>
#include <ranges>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

class queued_spinlock final {
public:
  class holder final {
    friend class queued_spinlock;

  public:
    explicit holder(queued_spinlock& host) : host_(host) {
      host_.acquire(this);
    }

    holder(const holder&) = delete;
    holder& operator=(const holder&) = delete;
    holder(holder&&) = delete;
    holder& operator=(holder&&) = delete;

    ~holder() { host_.release(this); }

  private:
    queued_spinlock& host_;
    std::atomic<holder*> next_ = nullptr;
    std::atomic<bool> holder_ = false;
  };

private:
  void acquire(holder* guard) noexcept {
    holder* prev = tail_.exchange(guard, std::memory_order::acq_rel);

    if (prev != nullptr) {
      prev->next_.store(guard, std::memory_order::release);

      while (!guard->holder_.load(std::memory_order::acquire)) {
        while (!guard->holder_.load(std::memory_order::relaxed))
          ;
      }
    } else {
      guard->holder_.store(true, std::memory_order::relaxed);
    }
  }

  void release(holder* guard) noexcept {
    if (guard->next_.load(std::memory_order::acquire) == nullptr) {
      holder* expected = guard;

      if (tail_.compare_exchange_strong(expected, nullptr,
                                        std::memory_order::release,
                                        std::memory_order::relaxed)) {
        return;
      }

      while (guard->next_.load(std::memory_order::acquire) == nullptr) {
        while (guard->next_.load(std::memory_order::relaxed) == nullptr)
          ;
      }
    }

    guard->next_.load(std::memory_order::relaxed)
        ->holder_.store(true, std::memory_order::release);
  }

private:
  std::atomic<holder*> tail_ = nullptr;
};

[[nodiscard]] inline auto next_file(fs::directory_iterator& it,
                                    fs::path& path) noexcept -> bool {
  static queued_spinlock lock;
  queued_spinlock::holder spin{lock};
  
  if (it == fs::directory_iterator{}) [[unlikely]] {
    return false;
  }

  path = it->path();
  ++it;
  return true;
}

template <typename T>
class el_vector {
  std::vector<T> vec_;
  std::size_t idx_ = 0;

public:
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;

  auto size() const noexcept -> std::size_t { return idx_; }

  auto begin() const noexcept -> const_iterator { return vec_.begin(); }
  auto end() const noexcept -> const_iterator { return vec_.begin() + idx_; }

  auto begin() noexcept -> iterator { return vec_.begin(); }
  auto end() noexcept -> iterator { return vec_.begin() + idx_; }

  auto vec() noexcept -> std::vector<T>& { return vec_; }

  auto operator[](std::size_t i) const noexcept -> T {
    return vec_[i];
  }

  void add(T val) {
    if (idx_ < vec_.size()) {
      vec_[idx_] = val;
    } else {
      vec_.push_back(val);
    }
    ++idx_;
  }

  void set_size(std::size_t n) noexcept { idx_ = n; }
};

static inline void parse_file_frugal(
    const fs::path& file,
    std::size_t k,
    std::size_t s,
    std::size_t t,
    el_vector<kmer_t>& dst
) {
  parser skp{file};
  const kmer_t mask = (1ULL << (2 * k)) - 1;

  for (const auto& seq : skp.sequences()) {
    kmer_t kmer = 0;
    kmer_t revc = 0;
    std::size_t base_c = 0;

    for (const char c : std::get<1>(seq)) {
      const auto base = base2(c);

      if (base < 0) {
        base_c = 0;
        continue;
      }

      kmer = (kmer << 2) | base;
      revc = (revc >> 2) | (static_cast<kmer_t>(3 - base) << (2 * (k - 1)));
      
      if (++base_c >= k) {
        const auto canon = std::min(kmer & mask, revc & mask);

        if (sync(canon, k, s, t)) {
          dst.add(canon);
        }
      }
    }
  }
}

class sequence_header final {
public:
  enum flags : std::uint8_t {};

private:
  static constexpr std::uint64_t lg_max_h = 32;
  static constexpr std::uint64_t len_bits = 12;
  static constexpr std::uint64_t len_offs = len_bits + lg_max_h;
  static constexpr std::uint64_t flg_offs = len_bits + len_offs;

  static constexpr std::uint64_t off_mask = (std::uint64_t{1} << len_offs) - 1;
  static constexpr std::uint64_t len_mask = ((std::uint64_t{1} << len_bits) - 1) << len_offs;
  static constexpr std::uint64_t flg_mask = ~(off_mask | len_mask);

  std::uint64_t data_;

  inline auto offset() const noexcept -> std::size_t {
    return static_cast<std::size_t>(data_ & off_mask);
  }

  inline auto length() const noexcept -> std::size_t {
    return static_cast<std::size_t>((data_ & len_mask) >> len_offs);
  }

public:
  sequence_header() = default;

  sequence_header(sequence_header&&) noexcept = default;
  sequence_header& operator=(sequence_header&&) noexcept = default;

  sequence_header(const sequence_header&) = delete;
  sequence_header& operator=(const sequence_header&) = delete;

  sequence_header(std::string& s, std::string&& l, flags f = {}) noexcept
      : data_(s.size())
  {
    data_ |= static_cast<std::uint64_t>(l.size()) << len_offs;
    data_ |= static_cast<std::uint64_t>(f) << flg_offs;
    s.append(std::move(l));
  }

  [[nodiscard]] auto name(const std::string& s) const noexcept -> std::string_view {
    return std::string_view{s}.substr(offset(), length());
  }

  [[nodiscard]] auto meta() const noexcept -> flags {
    return static_cast<flags>((data_ & flg_mask) >> flg_offs);
  }

  template<class Archive> void serialize(Archive& ar) { ar(data_); }
};

struct seq_view {
  std::size_t card;
  sequence_header file;

  seq_view(std::size_t c, std::string& s, std::string&& l, sequence_header::flags f = {})
      : card(c), file(s, std::move(l), f) {}
};

struct merge_view {
  std::size_t size;
  std::string path;
};

struct intermediate {
  std::vector<merge_view> parts;
  std::string names;
  std::vector<seq_view> seq;
};

constexpr std::string_view temp_dir = "tmp";

class file_cart {
  static constexpr std::size_t max_len = 8;
  static constexpr std::size_t tid_len = 2;
  static constexpr std::string_view def_tid = "000000XX";

public:
  class bag {
    friend class file_cart;

  public:
    inline auto next_filename() -> std::string {
      ++fid_[idx_];
      if (idx_ < 0 || fid_[idx_] > '~') {
        idx_ = max_len - tid_len - 1;
      }
      return fid_;
    }

  private:
    explicit bag(std::string&& tid) : fid_(std::move(tid)) {}

    int idx_ = static_cast<int>(max_len - tid_len - 1);
    std::string fid_;
  };

  [[nodiscard]] inline auto get_bag() noexcept -> bag {
    std::uint32_t hi = 0;
    std::uint32_t lo = 0;
    {
      queued_spinlock::holder qsh{tid_lock_};
      hi = next_hi_;
      lo = next_lo_;
      if (++next_lo_ > '~') {
        next_lo_ = '0';
        ++next_hi_; 
      }
    }
    std::string tid{def_tid};
    tid[max_len - tid_len] = static_cast<char>(hi);
    tid[max_len - tid_len + 1] = static_cast<char>(lo);
    return bag(std::move(tid));
  }

private:
  queued_spinlock tid_lock_;
  std::uint32_t next_hi_ = '0';
  std::uint32_t next_lo_ = '0';
};

[[nodiscard]] static inline auto extract_non_owning(
    fs::directory_iterator& dir,
    fs::path& file,
    std::size_t k,
    std::size_t s,
    std::size_t t,
    el_vector<kmer_t>& vec,
    std::string& labels,
    std::vector<seq_view>& seq,
    queued_spinlock& seq_lock
) -> bool {
  const std::size_t prev_size = vec.size();

  if (!next_file(dir, file)) [[unlikely]] {
    return false;
  }
  
  parse_file_frugal(file, k, s, t, vec);

  const auto start = vec.begin() + prev_size;

  boost::sort::spreadsort::spreadsort(start, vec.end());

  auto q = std::ranges::unique(start, vec.end()).begin();

  const auto card = static_cast<std::size_t>(q - start);
  {
    queued_spinlock::holder spin{seq_lock};
    seq.emplace_back(card, labels, file.filename().string());
  }

  vec.set_size(q - vec.begin());

  return true;
}

[[nodiscard]] static inline auto buffered_extract(
    fs::directory_iterator& dir,
    fs::path& file,
    std::size_t k,
    std::size_t s,
    std::size_t t,
    el_vector<kmer_t>& vec,
    std::size_t max_elems,
    std::string& labels,
    std::vector<seq_view>& seq,
    queued_spinlock& seq_lock
) -> bool {
  vec.set_size(0);

  while (vec.size() < max_elems) {
    if (!extract_non_owning(dir, file, k, s, t, vec, labels, seq, seq_lock)) {
      return false;
    }
  }

  boost::sort::spreadsort::spreadsort(vec.begin(), vec.end());

  auto q = std::ranges::unique(vec.begin(), vec.end()).begin();

  vec.set_size(q - vec.begin());

  return true;
}

inline void worker(
    fs::directory_iterator& dir,
    std::size_t k,
    std::size_t s,
    std::size_t t,
    std::size_t worker_limit,
    intermediate& im,
    queued_spinlock& seqlock,
    queued_spinlock& batlock,
    file_cart& cart
) {
  el_vector<kmer_t> kmer_buf;
  fs::path input_file_buf;

  file_cart::bag gen = cart.get_bag();
  std::ofstream os;

  std::size_t max_elems = worker_limit / (2 * sizeof(kmer_t));
  kmer_buf.vec().reserve(max_elems);

  while (true) {
    auto it = buffered_extract(dir, input_file_buf, k, s, t, kmer_buf,
                               max_elems, im.names, im.seq, seqlock);
    if (!it && kmer_buf.size() == 0) {
      break;
    }
    
    std::string label = gen.next_filename();
    {
      queued_spinlock::holder qsh{batlock};
      im.parts.emplace_back(kmer_buf.size(), label);
    }
    
    fs::path file = fs::path(std::string(temp_dir)) / label;

    os.clear();
    os.open(file, std::ios::binary | std::ios::trunc);
    if (!os) {
      throw std::runtime_error{"failed to open temporary file"};
    }

    os.write(reinterpret_cast<const char*>(kmer_buf.vec().data()), kmer_buf.size() * sizeof(kmer_t));
    if (!os) {
      throw std::runtime_error{"failed to write temporary file"};
    }

    if (!it) {
      break;
    }
  }
}

auto min_sum_of_filename_sizes(std::size_t f) -> std::size_t {
  return (4 + static_cast<std::size_t>(ceil(log10(static_cast<double>(f))))) * f;
}

inline auto convert_all_files(
    fs::path dir,
    std::size_t k,
    std::size_t s,
    std::size_t t,
    std::size_t mem_limit
) -> intermediate {
  intermediate repr;
  queued_spinlock seqlock;
  queued_spinlock batlock;
  file_cart cart;

  fs::directory_iterator iter{dir};

  std::size_t files = 0;
  for (const auto& p : iter) {
    (void)p;
    ++files;
  }

  std::size_t def_workers = tbb::info::default_concurrency();

  fs::create_directory(temp_dir);

  repr.names.reserve(min_sum_of_filename_sizes(files));
  repr.seq.reserve(files);
  repr.parts.reserve(def_workers);

  iter = fs::directory_iterator{dir};
  std::size_t lim = mem_limit / def_workers;

  tbb::parallel_for(
    tbb::blocked_range<std::size_t>(0, def_workers),
    [&](const auto& rng) {
      (void)rng;
      worker(iter, k, s, t, lim, repr, seqlock, batlock, cart);
    }
  );

  return repr;
}