#include "mergesort.h"

int main(int argc, char** argv) {
  assert(argc > 1);

  const auto ap = fs::absolute(argv[1]);
  assert(fs::is_directory(ap));

  convert_all_files(ap, 31, 0, 0, 4ULL * (1ULL << 32));

  return 0;
}