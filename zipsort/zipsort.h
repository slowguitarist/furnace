// Zipped iterator using two random-access iterators.
// To be until P3179R9 is implemented.

#pragma once

#include <iterator>
#include <execution>

template<class T, class U>
struct zip_value {
  T first;
  U second;
};

template<class T, class U>
class zip_reference {
public:
  T& first;
  U& second;

  using value_type = zip_value<std::remove_cv_t<T>, std::remove_cv_t<U>>;

  zip_reference(T& x, U& y) : first(x), second(y) {}

  operator value_type() && {
    return {std::move(first), std::move(second)};
  }

  zip_reference& operator=(const zip_reference& rhs) {
    first = rhs.first;
    second = rhs.second;
    return *this;
  }

  zip_reference& operator=(zip_reference&& rhs) {
    first = std::move(rhs.first);
    second = std::move(rhs.second);
    return *this;
  }

  zip_reference& operator=(const value_type& rhs) {
    first = rhs.first;
    second = rhs.second;
    return *this;
  }

  zip_reference& operator=(value_type&& rhs) {
    first = std::move(rhs.first);
    second = std::move(rhs.second);
    return *this;
  }

  friend void swap(zip_reference a, zip_reference b)
      noexcept(noexcept(std::ranges::swap(a.first, b.first)) &&
               noexcept(std::ranges::swap(a.second, b.second))) {
    std::ranges::swap(a.first, b.first);
    std::ranges::swap(a.second, b.second);
  }
};

template<class It1, class It2>
class zip_iterator {
  static_assert(std::is_lvalue_reference_v<std::iter_reference_t<It1>>);
  static_assert(std::is_lvalue_reference_v<std::iter_reference_t<It2>>);

  It1 it1_{};
  It2 it2_{};

public:
  using reference = zip_reference<
      std::remove_reference_t<std::iter_reference_t<It1>>,
      std::remove_reference_t<std::iter_reference_t<It2>>>;
  using value_type = typename reference::value_type;
  using difference_type = std::common_type_t<
      std::iter_difference_t<It1>,
      std::iter_difference_t<It2>>;
  using pointer = void;
  using iterator_category = std::random_access_iterator_tag;
  using iterator_concept = std::random_access_iterator_tag;

  zip_iterator() = default;
  zip_iterator(It1 a, It2 b) : it1_(a), it2_(b) {}

  reference operator*() const { return {*it1_, *it2_}; }
  reference operator[](difference_type n) const { return {it1_[n], it2_[n]}; }

  zip_iterator& operator++() { ++it1_; ++it2_; return *this; }
  zip_iterator operator++(int) { auto old = *this; ++*this; return old; }
  zip_iterator& operator--() { --it1_; --it2_; return *this; }
  zip_iterator operator--(int) { auto old = *this; --*this; return old; }

  zip_iterator& operator+=(difference_type n)
    { it1_ += n; it2_ += n; return *this; }
  zip_iterator& operator-=(difference_type n)
    { it1_ -= n; it2_ -= n; return *this; }

  friend zip_iterator operator+(zip_iterator it, difference_type n)
    { return it += n; }
  friend zip_iterator operator+(difference_type n, zip_iterator it)
    { return it += n; }
  friend zip_iterator operator-(zip_iterator it, difference_type n)
    { return it -= n; }

  friend difference_type operator-(const zip_iterator& a, const zip_iterator& b)
    { return a.it1_ - b.it1_; }

  friend bool operator==(const zip_iterator& a, const zip_iterator& b)
    { return a.it1_ == b.it1_ && a.it2_ == b.it2_; }
  friend bool operator<(const zip_iterator& a, const zip_iterator& b)
    { return a.it1_ < b.it1_; }
  friend bool operator>(const zip_iterator& a, const zip_iterator& b)
    { return b < a; }
  friend bool operator<=(const zip_iterator& a, const zip_iterator& b)
    { return !(b < a); }
  friend bool operator>=(const zip_iterator& a, const zip_iterator& b)
    { return !(a < b); }
};

template <typename Ep, typename Ct1, typename Ct2, typename F>
void zipsort(Ep&& ep, Ct1& c1, Ct2& c2, F cmp) {
  auto f1 = std::begin(c1);
  auto l1  = std::end(c1);
  auto f2 = std::begin(c2);
  auto l2  = std::end(c2);

  if (l1 - f1 != l2 - f2) {
    throw std::invalid_argument("zipsort: containers must have equal lengths.");
  }

  using It = zip_iterator<decltype(f1), decltype(f2)>;
  std::sort(std::forward<Ep>(ep), It{f1, f2}, It{l1, l2}, std::move(cmp));
}