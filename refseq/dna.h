// Quick and handy FASTA syncmer extractor.

#pragma once

#include <cassert>
#include <filesystem>
#include <generator>
#include <fstream>
#include <ranges>
#include <string>
#include <vector>

#include <boost/sort/spreadsort/spreadsort.hpp>

#include <sdsl/sd_vector.hpp>

#define RM_SEQ_DIR_ARG()                                                  \
  do {                                                                    \
    for (int iii = 1; iii < argc - 1; ++iii) argv[iii] = argv[iii + 1];   \
    argv[--argc] = nullptr;                                               \
  } while (0);

template<class T> void dust(const T&) {}

#define FORCE_LINK(x) { void* volatile ddd = &(x); dust(ddd); }

namespace fs = std::filesystem;

using kmer_t = std::uint64_t;
using kset_t = sdsl::sd_vector<>;

class parser {
public:
  using value_type = std::pair<std::string, std::string>;

  parser(const fs::path& path) : is_{path} {}

  auto sequences() -> std::generator<value_type> {
    if (!is_) {
      co_return;
    }
    
    std::string buf;
    value_type seq;

    while (std::getline(is_, buf) && (buf.empty() || buf[0] != '>'))
      ;  

    auto& name = std::get<0>(seq);
    auto& read = std::get<1>(seq);

    while (is_) {
      auto l = buf.size();

      if (buf[l - 1] == '\r') { l--; }

      name.assign(buf.begin() + 1, buf.begin() + l);
      read.clear();

      do {
        buf.clear();
        std::getline(is_, buf);
        if (buf.empty()) {
          continue;
        }

        if ((buf[0] != ';') && (buf[0] != '>')) {
          auto l = buf.size() - 1;
          auto end = (buf[l] == '\r') ? l : buf.size();
          read.insert(read.end(), buf.begin(), buf.begin() + end);
        } else if (buf[0] == '>') {
          co_yield std::move(seq);
          break;
        }
      } while (is_);

      if (!read.empty()) {
        co_yield std::move(seq);
      }
    }
  }

private:
  parser(const parser&) = delete;
  void operator=(const parser&) = delete;

  std::ifstream is_;
};

[[nodiscard]] static inline auto base2(char c) noexcept -> int {
  switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return -1;
  }
}

[[nodiscard]] static inline auto sync(kmer_t kmer,
                                      std::size_t k,
                                      std::size_t s,
                                      std::size_t t) noexcept -> bool {
  if (s == 0 || s >= k) {
    return true;
  }

  kmer_t smer_mask = (1U << (2 * s)) - 1;
  std::size_t num_smers = k - s + 1;

  std::size_t tmer_shift = 2 * (k - s - t);
  kmer_t tmer = (kmer >> tmer_shift) & smer_mask;

  for (std::size_t i = 0; i < num_smers; ++i) {
    kmer_t smer = kmer & smer_mask;
    if (smer < tmer) {
      return false;
    }
    kmer >>= 2;
  }

  return true;
}

static inline void parse_file(const fs::path& file,
                              std::size_t k,
                              std::size_t s,
                              std::size_t t,
                              std::vector<kmer_t>& dst) {
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
          dst.push_back(canon);
        }
      }
    }
  }
}

[[nodiscard]] static inline auto extract(const fs::path& file,
                                         std::size_t k,
                                         std::size_t s,
                                         std::size_t t) noexcept -> std::vector<kmer_t> {
  std::vector<kmer_t> vec;
  parser skp{file};

  std::error_code ec;
  auto hint = file_size(file, ec);

  if (ec) {
    vec.reserve(hint / sizeof(kmer_t));
  }

  parse_file(file, k, s, t, vec);

  return vec;
}

[[nodiscard]] static inline auto extract_sd(
    const fs::path& file,
    std::size_t k,
    std::size_t s,
    std::size_t t) noexcept -> std::pair<kset_t, std::size_t>
{
  std::vector<kmer_t> vec = extract(file, k, s, t);
  boost::sort::spreadsort::spreadsort(vec);

  auto uniq = std::ranges::unique(vec);
  const std::size_t card = uniq.begin() - vec.begin();

  return {kset_t{vec.begin(), uniq.begin()}, card};
}

[[nodiscard]] static inline auto assemblies(const char* dir) -> std::generator<fs::path> {
  using gotcha_tysm = fs::recursive_directory_iterator;

  const auto ap = fs::absolute(dir);
  assert(fs::is_directory(ap));

  for (const auto& f : gotcha_tysm(ap)) {
    co_yield f.path();
  }
}

using unique_range = std::ranges::borrowed_subrange_t<std::vector<kmer_t>&>;
using sorted_seq = std::tuple<std::vector<kmer_t>, unique_range, std::size_t>;

[[nodiscard]] static inline auto extract_master(const char* dir,
                                                std::size_t k,
                                                std::size_t s,
                                                std::size_t t) noexcept -> sorted_seq {
  std::vector<kmer_t> master;

  for (const auto path : assemblies(dir)) {
    parse_file(path, k, s, t, master);
  }

  boost::sort::spreadsort::spreadsort(master);

  auto uniq = std::ranges::unique(master);
  const std::size_t card = uniq.begin() - master.begin();

  return {std::move(master), uniq, card};
}

[[nodiscard]] static inline auto extract_all(
    const char* dir,
    std::size_t k,
    std::size_t s,
    std::size_t t) noexcept -> std::generator<std::vector<kmer_t>>
{
  for (const auto path : assemblies(dir)) {
    std::vector vec = extract(path, k, s, t);
    boost::sort::spreadsort::spreadsort(vec);
    vec.erase(std::ranges::unique(vec).begin(), vec.end());
    co_yield vec;
  }
}