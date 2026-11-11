#pragma once

#include "refseq/dna.h"
#include "definitions.h"
#include "utilities.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <limits>
#include <optional>
#include <queue>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_for.h>

namespace fs = std::filesystem;

[[nodiscard]] static inline auto extract_non_owning(
    fs::directory_iterator& dir,
    fs::path& file,
    std::size_t k,
    std::size_t s,
    std::size_t t,
    lazyvec<kmer_t>& vec,
    std::string& labels,
    std::vector<sview>& seq,
    queued_spinlock& seq_lock
) -> bool {
  const std::size_t prev_size = vec.size();

  if (!next_file(dir, file)) [[unlikely]] {
    return false;
  }
  
  simple_parse(file, k, s, t, vec);

  const auto start = vec.begin() + prev_size;

  boost::sort::spreadsort::spreadsort(start, vec.end());

  auto q = std::ranges::unique(start, vec.end()).begin();

  const auto card = static_cast<std::size_t>(q - start);
  {
    queued_spinlock::holder spin{seq_lock};
    seq.emplace_back(card, labels, file.filename().string());
  }

  vec.resize(q - vec.begin());

  return true;
}

[[nodiscard]] static inline auto buffered_extract(
    fs::directory_iterator& dir,
    fs::path& file,
    std::size_t k,
    std::size_t s,
    std::size_t t,
    lazyvec<kmer_t>& vec,
    std::size_t max_elems,
    std::string& labels,
    std::vector<sview>& seq,
    queued_spinlock& seq_lock
) -> bool {
  vec.resize(0);

  while (vec.size() < max_elems) {
    if (!extract_non_owning(dir, file, k, s, t, vec, labels, seq, seq_lock)) {
      return false;
    }
  }

  boost::sort::spreadsort::spreadsort(vec.begin(), vec.end());

  auto q = std::ranges::unique(vec.begin(), vec.end()).begin();

  vec.resize(q - vec.begin());

  return true;
}

inline void worker(
    fs::directory_iterator& dir,
    std::size_t k,
    std::size_t s,
    std::size_t t,
    std::size_t worker_limit,
    mergesort_components& im,
    queued_spinlock& seqlock,
    queued_spinlock& batlock,
    file_cart& cart
) {
  lazyvec<kmer_t> kbuf;
  fs::path input_file_buf;

  file_cart::bag gen = cart.get_bag();
  std::ofstream os;

  std::size_t max_elems = worker_limit / (2 * sizeof(kmer_t));
  kbuf.vec().reserve(max_elems);

  while (true) {
    auto it = buffered_extract(dir, input_file_buf, k, s, t, kbuf,
                               max_elems, im.names, im.seq, seqlock);
    if (!it && kbuf.size() == 0) {
      break;
    }
    
    std::string label = gen.next_filename();
    {
      queued_spinlock::holder qsh{batlock};
      im.parts.emplace_back(kbuf.size(), label);
    }
    
    fs::path file = get_temp_dir() / label;

    os.clear();
    os.open(file, std::ios::binary | std::ios::trunc);
    if (!os) {
      throw std::runtime_error{"failed to open temporary file"};
    }

    os.write(reinterpret_cast<const char*>(kbuf.vec().data()), kbuf.size() * sizeof(kmer_t));
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

[[nodiscard]] inline auto convert_all_files(
    fs::path dir,
    std::size_t k,
    std::size_t s,
    std::size_t t,
    std::size_t mem_limit,
    file_cart& cart
) -> mergesort_components {
  mergesort_components repr;
  queued_spinlock seqlock;
  queued_spinlock batlock;
  
  fs::directory_iterator iter{dir};

  std::size_t files = 0;
  for (const auto& p : iter) {
    (void)p;
    ++files;
  }

  const std::size_t def_workers = tbb::info::default_concurrency();
  const std::size_t work = std::min(def_workers, files);

  fs::create_directory(get_temp_dir());

  repr.names.reserve(min_sum_of_filename_sizes(files));
  repr.seq.reserve(files);
  repr.parts.reserve(work);

  iter = fs::directory_iterator{dir};
  std::size_t lim = mem_limit / def_workers;

  tbb::parallel_for(
    tbb::blocked_range<std::size_t>(0, work),
    [&](const auto& rng) {
      (void)rng;
      worker(iter, k, s, t, lim, repr, seqlock, batlock, cart);
    }
  );

  return repr;
}


struct mstream {
  static constexpr std::size_t buf_size = 512ULL * 1024ULL;
  static constexpr auto eof = std::numeric_limits<std::size_t>::max();

  std::ifstream ifs;
  std::span<kmer_t> buf;
  std::size_t idx = 0;
  std::size_t size = 0;

  void init(mview&& mv, std::span<kmer_t> bp) {
    assert(bp.size() == buf_size);
    buf = bp;
    ifs.open(mv.path, std::ios::binary);
    if (!ifs) {
      throw std::runtime_error{"could not open input merge " + mv.path};
    }
    refill();
  }

  void refill() {
    if (idx == eof) {
      return;
    }
    ifs.read(reinterpret_cast<char*>(buf.data()), buf.size() * sizeof(kmer_t));
    size = ifs.gcount() / sizeof(kmer_t);
    if (size == 0) {
      idx = eof;
      ifs.close();
    } else {
      idx = 0;
    }
  }

  [[nodiscard]] auto next() -> std::optional<kmer_t> {
    if (idx == eof) {
      return std::nullopt;
    }
    kmer_t res = buf[idx];
    if (++idx >= size) {
      refill();
    }
    return idx != eof ? std::make_optional(res) : std::nullopt; 
  }

  explicit operator bool() {
    return idx != eof;
  }
};

inline void k_way_loser_merge_non_owning(
    std::span<mview> parts,
    std::vector<kmer_t>& sbuf,
    std::vector<kmer_t>& outbuf,
    std::vector<mstream>& mstr,
    std::priority_queue<loser>& minh,
    file_cart::bag& gen
) {
  std::size_t gidx = 0;
  std::size_t osize = 0;

  const std::string spell = gen.next_filename();
  std::ofstream ofs(spell, std::ios::binary | std::ios::trunc);
  if (!ofs) {
    throw std::runtime_error{"could not open output merge " + spell};
  }

  for (auto [part, stream] : std::views::zip(parts, mstr)) {
    stream.init(std::move(part), std::span{sbuf.begin() + gidx, mstream::buf_size});
    gidx += mstream::buf_size;
  }

  gidx = 0;

  auto flush = [&]() {
    osize += gidx;
    ofs.write(reinterpret_cast<const char*>(outbuf.data()), gidx * sizeof(kmer_t));
    gidx = 0;
  };

  auto fetch = [&](std::size_t idx) {
    auto opt = mstr[idx].next();
    if (opt.has_value()) {
      minh.emplace(*opt, idx);
    }
  };

  for (std::size_t i = 0; i < mstr.size(); ++i) {
    fetch(i);
  }

  while (!minh.empty()) {
    const kmer_t dog = minh.top().val;
    do {
      std::size_t idx = minh.top().idx;
      minh.pop();
      fetch(idx);
    } while (!minh.empty() && minh.top().val == dog);

    outbuf[gidx] = dog;
    if (++gidx >= outbuf.size()) {
      flush();
    }
  }

  if (gidx > 0) {
    flush();
  }
  if (!ofs) {
    throw std::runtime_error{"failed writing local merge"};
  }

  parts.front() = {osize, std::move(spell)};
}