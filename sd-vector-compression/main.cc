#include "refseq/dna.h"

#include <cstddef>
#include <print>
#include <sstream>

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

int main(int argc, char** argv) {
  assert(argc > 1);

  std::size_t sum_sv = 0;
  std::size_t sum_sd = 0;
  std::size_t with_dup = 0;

  for (auto&& seq : extract_all(argv[1], 31, 0, 0)) {
    with_dup += seq.size();

    kset_t sd{seq.begin(), seq.end()};
  
    std::ostringstream ssv;
    {
      cereal::BinaryOutputArchive svar{ssv};
      svar(seq);
    }

    std::ostringstream ssd;
    {
      cereal::BinaryOutputArchive sdar{ssd};
      sdar(sd);
    }

    const std::size_t svb = ssv.view().size();
    const std::size_t sdb = ssd.view().size();
    const double ratio = static_cast<double>(svb) / static_cast<double>(sdb);

    sum_sv += svb;
    sum_sd += sdb;

    std::println("vec: {}, sdvec: {}, ratio: {}", svb, sdb, ratio);
  }

  auto [vec, uniq, card] = extract_master(argv[1], 31, 0, 0);

  vec.erase(uniq.begin(), vec.end());

  kset_t msd{vec.begin(), vec.end()};

  std::ostringstream smv;
  {
    cereal::BinaryOutputArchive smvar{smv};
    smvar(vec);
  }

  std::ostringstream smd;
  {
    cereal::BinaryOutputArchive smdar{smd};
    smdar(msd);
  }

  const std::size_t smvb = smv.view().size();
  const std::size_t smdb = smd.view().size();
  const double svr = static_cast<double>(sum_sv) / static_cast<double>(smvb);
  const double sdr = static_cast<double>(sum_sd) / static_cast<double>(smdb);

  std::println("----");
  std::println("dup: {}, mas: {}, rem: {}", with_dup, card, with_dup - card);
  std::println("vdup: {}, vmas: {}, ratio: {}", sum_sv, smvb, svr);
  std::println("ddup: {}, dmas: {}, ratio: {}", sum_sd, smdb, sdr);
}