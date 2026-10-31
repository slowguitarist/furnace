#include "refseq/dna.h"

#include <filesystem>
#include <print>

#include <benchmark/registration.h>
#include <benchmark/state.h>

#include <boost/sort/spreadsort/spreadsort.hpp>

namespace fs = std::filesystem;

using rdi = fs::recursive_directory_iterator;

int main(int argc, char** argv) {
  benchmark::MaybeReenterWithoutASLR(argc, argv);

  if (argc < 2) {
    std::println("argc {}", argc);
    return 1;
  }

  const fs::path f = fs::absolute(argv[1]);

  // Note to self: pass path as "$PWD/refseq/fasta"
  if (!fs::is_directory(f)) {
    std::println("Not dir: {}", argv[1]);
    return 1337;
  }

  for (int i = 1; i < argc - 1; ++i) argv[i] = argv[i + 1];
  argv[--argc] = nullptr;

  static constexpr std::size_t k = 31;
  static constexpr std::size_t s = 0;
  static constexpr std::size_t t = 0;

  for (const auto& file : rdi(f)) {
    const auto path = file.path();

    if (!file.is_regular_file()) {
      std::println("Not file: {}", path);
      return 6174;
    }

    std::vector<kmer_t> data = extract(path, k, s, t);
    std::vector<kmer_t> copy = data;

    const std::string radix = path.string() + "-boost";
    const std::string stand = path.string() + "-std";

    benchmark::RegisterBenchmark(radix.c_str(), [data](benchmark::State& st) mutable {
      for (auto _ : st) {
        boost::sort::spreadsort::spreadsort(data);
      }
    });

    benchmark::RegisterBenchmark(stand.c_str(), [copy](benchmark::State& st) mutable {
      for (auto _ : st) {
        std::ranges::sort(copy);
      }
    });
  }

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
}