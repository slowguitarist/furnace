#include "refseq/dna.h"

#include <cassert>
#include <optional>
#include <ranges>
#include <vector>

#include <bbhash/bbhash.h>

#include <benchmark/registration.h>
#include <benchmark/state.h>

auto test_seq(std::vector<kmer_t> master,
              unique_range uniq,
              std::size_t card,
              std::size_t size) -> std::pair<std::vector<kmer_t>, std::vector<kmer_t>>
{
  assert(size > 0);
  assert(card >= size * 2);

  std::size_t d = (card - 1) / size;

  std::vector<kmer_t> present;
  std::vector<kmer_t> absent;

  present.reserve(size);
  absent.reserve(size);

  for (auto j = master.begin() + 1; j < uniq.begin(); j += d) {
    present.push_back(*j);
    absent.push_back((*j + *(j - 1)) / 2);
  }

  return {std::move(present), std::move(absent)};
}

int main(int argc, char** argv) {
  benchmark::MaybeReenterWithoutASLR(argc, argv);
  assert(argc > 1);

  auto [vec, uniq, card] = extract_master(argv[1], 31, 0, 0);

  kset_t kset{vec.begin(), uniq.begin()};

  kset_t::rank_1_type rank(&kset);
  kset_t::select_1_type select(&kset);

  bbh::bbhash<kmer_t> hash(std::ranges::subrange(vec.begin(), uniq.begin()));

  auto [pr, ab] = test_seq(vec, uniq, card, 1337);

  benchmark::RegisterBenchmark("rank-present", [&](benchmark::State& st) {
    for (auto _ : st) {
      for (const kmer_t kmer : pr) {
        const std::size_t r = rank(kmer);
        auto op = (r < card && select(r + 1) == kmer) ? std::make_optional(r) : std::nullopt;
        assert(op.has_value());
      }
    }
  });

  benchmark::RegisterBenchmark("hash-present", [&](benchmark::State& st) {
    for (auto _ : st) {
      for (const kmer_t kmer : pr) {
        if (kmer < card && kset[kmer] == 1) {
          auto op = hash.find(kmer);
          assert(op.has_value());
        }
      } 
    }
  });

  benchmark::RegisterBenchmark("rank-maybe-absent", [&](benchmark::State& st) {
    for (auto _ : st) {
      for (const kmer_t kmer : ab) {
        const std::size_t r = rank(kmer);
        auto op = (r < card && select(r + 1) == kmer) ? std::make_optional(r) : std::nullopt;
        FORCE_LINK(op);
      }
    }
  });

  benchmark::RegisterBenchmark("hash-maybe-absent", [&](benchmark::State& st) {
    for (auto _ : st) {
      for (const kmer_t kmer : ab) {
        if (kmer < card && kset[kmer] == 1) {
          auto op = hash.find(kmer);
          FORCE_LINK(op);
        }
      } 
    }
  });

  RM_SEQ_DIR_ARG();

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();

  return 0;
}