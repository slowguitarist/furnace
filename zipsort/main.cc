#include "zipsort.h"

#include <execution>
#include <memory>
#include <print>
#include <vector>

struct Copy {
  std::uint64_t w;
  double x;
  int p;
};

template <>
struct std::formatter<Copy> : std::formatter<std::string_view> {
  auto format(const Copy& c, auto& ctx) const {
    return std::format_to(ctx.out(), "{{w:{}, x:{:.2f}, p:{}}}", c.w, c.x, c.p);
  }
};

struct Move {
  std::int8_t y;
  std::unique_ptr<std::size_t> z;
};

template <>
struct std::formatter<Move> : std::formatter<std::string_view> {
  auto format(const Move& m, auto& ctx) const {
    if (m.z) {
      return std::format_to(ctx.out(), "{{y:{}, z:{}}}", static_cast<int>(m.y), *m.z);
    } else {
      return std::format_to(ctx.out(), "{{y:{}, z:null}}", static_cast<int>(m.y));
    }
  }
};

int main() {
  std::vector<int>    v1 = {5,   7,   3,   1,   8,   0,   0};
  std::vector<double> v2 = {4.3, 1.5, 0.2, 7.8, 1.1, 2.3, 7.8};

  zipsort(std::execution::par, v1, v2, [](const auto& a, const auto& b) {
    return a.second > b.second;
  });

  std::println("{}", v1);
  std::println("{}", v2);

  constexpr std::size_t N = 1'000'000;

  std::vector<Copy> v3;
  std::vector<Move> v4;

  v3.reserve(N);
  v4.reserve(N);

  for (std::size_t i = 0; i < N; ++i) {
    v3.emplace_back(i, static_cast<double>(i & 2047), -static_cast<int>(i));
    v4.emplace_back(static_cast<std::int8_t>(i & 127), std::make_unique<std::size_t>(i * 2));
  }

  std::println("Init:");
  std::println("{}, {}, {}, {}", v3[0], v3[1], v3[2], v3[3]);
  std::println("{}, {}, {}, {}", v4[0], v4[1], v4[2], v4[3]);

  zipsort(std::execution::par, v3, v4, [](const auto& a, const auto& b) {
    return a.first.x > b.first.x;
  });

  std::println("Sort:");
  std::println("{}, {}, {}, {}", v3[0], v3[1], v3[2], v3[3]);
  std::println("{}, {}, {}, {}", v4[0], v4[1], v4[2], v4[3]);

  return 0;
}