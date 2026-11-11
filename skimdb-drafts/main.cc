#include <cassert>
#include <filesystem>

namespace fs = std::filesystem;

int main(int argc, char** argv) {
  assert(argc > 1);

  const fs::path ap = fs::absolute(argv[1]);
  assert(fs::is_directory(ap));

  // ...

  return 0;
}