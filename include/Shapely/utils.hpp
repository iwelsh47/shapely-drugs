#pragma once

#include <vector>

namespace shapely {
template <typename FloatT>
FloatT accurate_sum(const std::vector<FloatT>& values) {
  // Implements the msum function from http://code.activestate.com/recipes/393090/
  std::vector<FloatT> partials;
  partials.reserve(1 + values.size() / 10);
  for (uint64_t i = 0; i < values.size(); ++i) {
    FloatT x = values[i];
    uint64_t j = 0;
    for (FloatT y : partials) {
      if (std::abs(x) < std::abs(y)) { std::swap(x, y); }
      FloatT high = x + y;
      FloatT low = y - (high - x);
      if (low != 0) {
        partials[j++] = low;
      }
      x = high;
    }
    partials.resize(++j);
    partials.back() = x;
  }
  FloatT sum = 0;
  for (const FloatT& v : partials) { sum += v; }
  return sum;
}

template <int64_t N, typename T>
inline constexpr T pow_n(const T& x) {
  if constexpr(N == 0) { return 1; }
  if constexpr(N < 0) { return 1 / pow_n<-N>(x); }
  if constexpr(N & 1) { return x * pow_n<N - 1>(x); }
  else {
    T x2 = pow_n<N / 2>(x);
    return x2 * x2;
  }
}

}
