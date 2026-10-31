#include "refseq/dna.h"

#include <cassert>

int main(int argc, char** argv) {
  assert(argc > 1);

  for (auto f : assemblies(argv[1])) {
    auto [bv, card] = extract_sd(f, 31, 0, 0);
    assert(card > 0);
  }

  return 0;
}