#pragma once

#include "refseq/dna.h"

#include <algorithm>
#include <cstddef>
#include <fcntl.h>
#include <filesystem>
#include <ranges>
#include <vector>

class queued_spinlock final {
public:
  class holder final {
    friend class queued_spinlock;

  public:
    explicit holder(queued_spinlock& host) : host_(host) {
      host_.acquire(this);
    }

    holder(const holder&) = delete;
    holder& operator=(const holder&) = delete;
    holder(holder&&) = delete;
    holder& operator=(holder&&) = delete;

    ~holder() { host_.release(this); }

  private:
    queued_spinlock& host_;
    std::atomic<holder*> next_ = nullptr;
    std::atomic<bool> holder_ = false;
  };

private:
  void acquire(holder* guard) noexcept {
    holder* prev = tail_.exchange(guard, std::memory_order::acq_rel);

    if (prev != nullptr) {
      prev->next_.store(guard, std::memory_order::release);

      // todo cpu relax
      while (!guard->holder_.load(std::memory_order::acquire)) {
        while (!guard->holder_.load(std::memory_order::relaxed))
          ;
      }
    } else {
      guard->holder_.store(true, std::memory_order::relaxed);
    }
  }

  void release(holder* guard) noexcept {
    if (guard->next_.load(std::memory_order::acquire) == nullptr) {
      holder* expected = guard;

      if (tail_.compare_exchange_strong(expected, nullptr,
                                        std::memory_order::release,
                                        std::memory_order::relaxed)) {
        return;
      }

      while (guard->next_.load(std::memory_order::acquire) == nullptr) {
        while (guard->next_.load(std::memory_order::relaxed) == nullptr)
          ;
      }
    }

    guard->next_.load(std::memory_order::relaxed)
        ->holder_.store(true, std::memory_order::release);
  }

private:
  std::atomic<holder*> tail_ = nullptr;
};


[[nodiscard]] inline auto next_file(fs::directory_iterator& it,
                                    fs::path& path) noexcept -> bool {
  static queued_spinlock lock;

  queued_spinlock::holder spin{lock};
  
  if (it == fs::directory_iterator{}) [[unlikely]] {
    return false;
  }

  path = it->path();
  ++it;
  
  return true;
}

template <typename T>
class frugal_vector {
  std::vector<T> vec_;
  std::size_t idx_ = 0;

public:
  using iterator = std::vector<T>::iterator;
  using const_iterator = std::vector<T>::const_iterator;

  auto size() const noexcept -> std::size_t { return idx_; }

  auto begin() const noexcept -> const_iterator { return vec_.begin(); }
  auto end() const noexcept -> const_iterator { return vec_.begin() + idx_; }

  auto begin() noexcept -> iterator { return vec_.begin(); }
  auto end() noexcept -> iterator { return vec_.begin() + idx_; }

  auto vec() const noexcept -> std::vector<T>& { return vec_; }

  auto operator[](std::size_t i) const noexcept -> T {
    return vec_[i];
  }

  void add(T val) {
    if (idx_ < vec_.size()) {
      vec_[idx_] = val;
    } else {
      vec_.push_back(val);
    }
    ++idx_;
  }

  void set_size(std::size_t n) noexcept { idx_ = n; }
};

static inline void parse_file_frugal(
    const fs::path& file,
    std::size_t k,
    std::size_t s,
    std::size_t t,
    frugal_vector<kmer_t>& dst) {
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
          dst.add(canon);
        }
      }
    }
  }
}

[[nodiscard]] static inline auto extract_non_owning(
    fs::directory_iterator& dir,
    fs::path& file,
    std::size_t k,
    std::size_t s,
    std::size_t t,
    frugal_vector<kmer_t>& vec
) -> bool {
  const std::size_t prev_size = vec.size();

  if (!next_file(dir, file)) [[unlikely]] {
    return false;
  }
  
  parse_file_frugal(file, k, s, t, vec);

  const auto start = vec.begin() + prev_size;

  boost::sort::spreadsort::spreadsort(start, vec.end());

  auto q = std::ranges::unique(start, vec.end()).begin();

  const auto card = static_cast<std::size_t>(q - start);
  // push card, path

  vec.set_size(q - vec.begin());

  return true;
}

[[nodiscard]] static inline auto buffered_extract(
    fs::directory_iterator& dir,
    fs::path& file,
    std::size_t k,
    std::size_t s,
    std::size_t t,
    frugal_vector<kmer_t>& vec,
    std::size_t worker_limit
) -> frugal_vector<kmer_t>::iterator {
  vec.set_size(0);

  // Conservatively half the quota to account for additional space used by
  // spreadsort
  while (vec.size() < worker_limit / 2) {
    if (!extract_non_owning(dir, file, k, s, t, vec)) {
      break;
    }
  }

  boost::sort::spreadsort::spreadsort(vec.begin(), vec.end());

  return std::ranges::unique(vec.begin(), vec.end()).begin();
}