#pragma once

#include "refseq/dna.h"

#include <cstdint>
#include <string>
#include <type_traits>
#include <vector>

#include <bbhash/bbhash.h>

template <typename T>
  requires std::is_trivially_destructible_v<T>
class lazyvec {
  std::vector<T> vec_;
  std::size_t idx_ = 0;

public:
  using value_type = T;
  using size_type = std::size_t;
  using iterator = typename std::vector<T>::iterator;
  using const_iterator = typename std::vector<T>::const_iterator;
  using reference = T&;
  using const_reference = const T&;

  lazyvec() = default;
  lazyvec(std::size_t size) : vec_(size) {}

  size_type size() const noexcept { return idx_; }
  bool empty() const noexcept { return idx_ == 0; }

  iterator begin() noexcept { return vec_.begin(); }
  iterator end() noexcept { return vec_.begin() + idx_; }

  const_iterator begin() const noexcept { return vec_.begin(); }
  const_iterator end() const noexcept { return vec_.begin() + idx_; }

  reference front() { return vec_[0]; }
  const_reference front() const { return vec_[0]; }

  reference back() { return vec_[idx_]; }
  const_reference back() const { return vec_[idx_]; }

  reference operator[](size_type i) noexcept { return vec_[i]; }
  const_reference operator[](size_type i) const noexcept { return vec_[i]; }

  void push_back(const T& val) {
    if (idx_ < vec_.size()) {
      vec_[idx_] = val;
    } else {
      vec_.push_back(val);
    }
    ++idx_;
  }

  void push_back(T&& val) {
    if (idx_ < vec_.size()) {
      vec_[idx_] = std::move(val);
    } else {
      vec_.push_back(std::move(val));
    }
    ++idx_;
  }

  template <typename... Args>
  reference emplace_back(Args&&... args) {
    if (idx_ < vec_.size()) {
      vec_[idx_] = T{std::forward<Args>(args)...};
    } else {
      vec_.emplace_back(std::forward<Args>(args)...);
    }
    return vec_[idx_++];
  }

  void pop_back() noexcept {
    --idx_;
  }

  void resize(size_type n) {
    idx_ = n;
    if (n > vec_.size()) {
      vec_.resize(n);
    }
  }

  std::vector<T>& vec() noexcept { return vec_; }
};

class label final {
public:
  enum flags : std::uint8_t {};

private:
  static constexpr std::uint64_t lg_max_h = 32;
  static constexpr std::uint64_t len_bits = 12;
  static constexpr std::uint64_t len_offs = len_bits + lg_max_h;
  static constexpr std::uint64_t flg_offs = len_bits + len_offs;

  static constexpr std::uint64_t off_mask = (std::uint64_t{1} << len_offs) - 1;
  static constexpr std::uint64_t len_mask = ((std::uint64_t{1} << len_bits) - 1) << len_offs;
  static constexpr std::uint64_t flg_mask = ~(off_mask | len_mask);

  std::uint64_t data_;

  inline auto offset() const noexcept -> std::size_t {
    return static_cast<std::size_t>(data_ & off_mask);
  }

  inline auto length() const noexcept -> std::size_t {
    return static_cast<std::size_t>((data_ & len_mask) >> len_offs);
  }

public:
  label() = default;

  label(label&&) noexcept = default;
  label& operator=(label&&) noexcept = default;

  label(const label&) = delete;
  label& operator=(const label&) = delete;

  label(std::string& s, std::string&& l, flags f = {}) noexcept
      : data_(s.size())
  {
    data_ |= static_cast<std::uint64_t>(l.size()) << len_offs;
    data_ |= static_cast<std::uint64_t>(f) << flg_offs;
    s.append(std::move(l));
  }

  [[nodiscard]] auto name(const std::string& s) const noexcept -> std::string_view {
    return std::string_view{s}.substr(offset(), length());
  }

  [[nodiscard]] auto meta() const noexcept -> flags {
    return static_cast<flags>((data_ & flg_mask) >> flg_offs);
  }

  template<class Archive> void serialize(Archive& ar) { ar(data_); }
};

class alignas(16) kview final {
private:
  std::atomic<std::uint64_t> offset_flags_ = 0;

public:
  enum io_status : std::uint64_t {
    External = std::uint64_t{0},
    Pending  = std::uint64_t{1} << 62,
    Present  = std::uint64_t{1} << 63,
  };

  static constexpr std::uint64_t state_mask =
      static_cast<std::uint64_t>(Pending) |
      static_cast<std::uint64_t>(Present);

  std::uint32_t size = 0;
  std::uint32_t card = 0;

  kview() = default;

  kview(std::uint64_t off, std::uint64_t sz, std::uint64_t cd) noexcept
      : offset_flags_(off),
        size(static_cast<std::uint32_t>(sz)),
        card(static_cast<std::uint32_t>(cd)) {}

  kview(kview&& other) noexcept
      : offset_flags_(other.offset_flags_.load(std::memory_order::relaxed)),
        size(other.size),
        card(other.card) {}

  kview& operator=(kview&& other) noexcept {
    if (this != &other) {
      offset_flags_.store(
          other.offset_flags_.load(std::memory_order::relaxed),
          std::memory_order::relaxed);
      size = other.size;
      card = other.card;
    }
    return *this;
  }

  kview(const kview&) = delete;
  kview& operator=(const kview&) = delete;

  static auto io_to_str(io_status st) noexcept -> std::string_view {
    switch (st) {
      case External: return "External";
      case Pending:  return "Pending";
      case Present:  return "Present";
      default:       return "Unknown";
    }
  } 

  auto offset() const noexcept -> std::size_t {
    return static_cast<std::size_t>(
        offset_flags_.load(std::memory_order::relaxed) & ~state_mask);
  }

  void shift_offset(std::uint64_t amt) noexcept {
    offset_flags_ = std::min(Pending - 1, offset_flags_ + amt);
  }

  auto state() const noexcept -> io_status {
    return static_cast<io_status>(
        offset_flags_.load(std::memory_order_acquire) & state_mask);
  }

  void set_state(io_status st) noexcept {
    std::uint64_t old = offset_flags_.load(std::memory_order::relaxed);
    std::uint64_t des = (old & ~state_mask) | static_cast<std::uint64_t>(st);

    while (!offset_flags_.compare_exchange_weak(
        old, des,
        std::memory_order::release,
        std::memory_order::relaxed)) {
      des = (old & ~state_mask) | static_cast<std::uint64_t>(st);
    }
  }
};

/* For mergesort */

struct sview {
  std::size_t card;
  label file;

  sview(std::size_t c, std::string& s, std::string&& l, label::flags f = {})
      : card(c), file(s, std::move(l), f) {}
};

struct mview {
  std::size_t size;
  std::string path;

  auto operator<(const mview& k) const noexcept -> bool {
    return size < k.size;
  }
};

struct mergesort_components {
  std::vector<mview> parts;
  std::string names;
  std::vector<sview> seq;
};

struct loser {
  kmer_t val;
  std::size_t idx;

  auto operator<(const loser& k) const noexcept -> bool {
    return val > k.val;
  }
};

/* For external storage */

struct indexed_stream final {
  std::ostringstream str;
  std::vector<std::reference_wrapper<kview>> idx;

  indexed_stream() : str(std::ios::out | std::ios::binary), idx() {}

  ~indexed_stream() = default;
};

// Fake "database" representation for prototyping purposes
struct fakedb {
  std::size_t unique;
  kmer_t maximum;
  std::vector<label> sequences;
  std::string labels;
  std::vector<std::vector<std::uint32_t>> data;
  bbh::bbhash<kmer_t> hash;
};