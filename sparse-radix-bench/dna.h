// Quick and handy FASTA syncmer extractor.
// Does not sort / deduplicate because this is being benchmarked.

#pragma once

#include <filesystem>
#include <generator>
#include <fstream>
#include <string>
#include <vector>

using kmer_t = std::uint64_t;

class parser {
public:
  using value_type = std::pair<std::string, std::string>;

  parser(const std::filesystem::path& path) : is_{path} {}

  auto sequences() -> std::generator<value_type> {
    if (!is_) { co_return; }
    
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
        if (buf.empty()) { continue; } 

        if ((buf[0] != ';') && (buf[0] != '>')) {
          auto l = buf.size() - 1;
          auto end = (buf[l] == '\r') ? l : buf.size();
          read.insert(read.end(), buf.begin(), buf.begin() + end);
        } else if (buf[0] == '>') {
          co_yield std::move(seq);
          break;
        }
      } while (is_);

      if (!read.empty()) { co_yield std::move(seq); }
    }
  }

private:
  parser(const parser&) = delete;
  void operator=(const parser&) = delete;

  std::ifstream is_;
};

[[nodiscard]] inline auto base2(char c) noexcept -> int {
  switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return -1;
  }
}

[[nodiscard]] inline auto sync(kmer_t kmer,
                               std::size_t k,
                               std::size_t s,
                               std::size_t t) noexcept -> bool {
  if (s == 0 || s >= k) { return true; }

  kmer_t smer_mask = (1U << (2 * s)) - 1;
  std::size_t num_smers = k - s + 1;

  std::size_t tmer_shift = 2 * (k - s - t);
  kmer_t tmer = (kmer >> tmer_shift) & smer_mask;

  for (std::size_t i = 0; i < num_smers; ++i) {
    kmer_t smer = kmer & smer_mask;
    if (smer < tmer) { return false; }
    kmer >>= 2;
  }

  return true;
}

[[nodiscard]] inline auto extract(const std::filesystem::path& file,
                                  std::size_t k,
                                  std::size_t s,
                                  std::size_t t) noexcept -> std::vector<kmer_t> {
  std::vector<kmer_t> vec;
  parser skp{file};

  std::error_code ec;
  auto hint = file_size(file, ec);

  if (ec) vec.reserve(hint / sizeof(kmer_t));

  const kmer_t mask = (1ULL << (2 * k)) - 1;

  for (const auto& seq : skp.sequences()) {
    kmer_t kmer = 0;
    kmer_t revc = 0;
    std::size_t base_c = 0;

    for (const char c : std::get<1>(seq)) {
      const auto base = base2(c);
      if (base < 0) { base_c = 0; continue; }

      kmer = (kmer << 2) | base;
      revc = (revc >> 2) | (static_cast<kmer_t>(3 - base) << (2 * (k - 1)));
      
      if (++base_c >= k) {
        const auto canon = std::min(kmer & mask, revc & mask);
        if (sync(canon, k, s, t)) { vec.push_back(canon); }
      }
    }
  }

  return vec;
}