#pragma once

#include "refseq/dna.h"
#include "definitions.h"

#include <atomic>
#include <emmintrin.h>
#include <semaphore>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;

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

      while (!guard->holder_.load(std::memory_order::acquire)) {
        _mm_pause();
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
        _mm_pause();
      }
    }

    guard->next_.load(std::memory_order::relaxed)
        ->holder_.store(true, std::memory_order::release);
  }

private:
  std::atomic<holder*> tail_ = nullptr;
};

class file_cart {
  static constexpr std::size_t max_len = 8;
  static constexpr std::size_t tid_len = 2;
  static constexpr std::string_view def_tid = "000000XX";

public:
  class bag {
    friend class file_cart;

  public:
    inline auto next_filename() -> std::string {
      ++fid_[idx_];
      if (--idx_ < 0 || fid_[idx_] > '~') {
        idx_ = max_len - tid_len - 1;
      }
      return fid_;
    }

  private:
    explicit bag(std::string&& tid) : fid_(std::move(tid)) {}

    int idx_ = static_cast<int>(max_len - tid_len - 1);
    std::string fid_;
  };

  [[nodiscard]] inline auto get_bag() noexcept -> bag {
    std::uint32_t hi = 0;
    std::uint32_t lo = 0;
    {
      queued_spinlock::holder qsh{tid_lock_};
      hi = next_hi_;
      lo = next_lo_;
      if (++next_lo_ > '~') {
        next_lo_ = '0';
        ++next_hi_; 
      }
    }
    std::string tid{def_tid};
    tid[max_len - tid_len] = static_cast<char>(hi);
    tid[max_len - tid_len + 1] = static_cast<char>(lo);
    return bag(std::move(tid));
  }

private:
  queued_spinlock tid_lock_;
  std::uint32_t next_hi_ = '0';
  std::uint32_t next_lo_ = '0';
};

class self_checkout final {
public:
  self_checkout(std::uint32_t count) 
      : tickets_(count),
        page_(count),
        stack_(std::make_unique<int[]>(count)) {
    for (std::uint32_t i = 0; i < count; ++i) {
      stack_[i] = i;
    }
  }

  int take() noexcept {
    // Can block if no tickets are available. "Help is on the way."
    tickets_.acquire();
    queued_spinlock::holder visit{kiosk_};
    return stack_[--page_];
  }

  void give(int ticket) noexcept {
    {
      queued_spinlock::holder visit{kiosk_};
      stack_[page_++] = ticket;
    }
    tickets_.release();
  }

private:
  queued_spinlock kiosk_;
  std::counting_semaphore<> tickets_;
  std::uint32_t page_;
  std::unique_ptr<int[]> stack_;
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

inline auto get_temp_dir() -> fs::path {
  if (const char* env_dir = std::getenv("BUILD_WORKING_DIRECTORY")) {
    return fs::path(env_dir) / "tmp";
  }
  return fs::path("tmp");
}

static inline void simple_parse(
    const fs::path& file,
    std::size_t k,
    std::size_t s,
    std::size_t t,
    lazyvec<kmer_t>& dst
) {
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
          dst.push_back(canon);
        }
      }
    }
  }
}