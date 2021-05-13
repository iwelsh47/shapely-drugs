#include <algorithm>
#include <cstdint>
#include <cmath>
#include <vector>
#include <Shapely/cubic_spline.hpp>

namespace shapely {

CubicSpline::CubicSpline(const std::vector<float>& x, const std::vector<float>& y)
: points(x), values(y), low(y.front()), high(y.back()) {
  // Perform snaity checks
  if (x.size() < 2) {
    throw std::range_error("CubicSpline requires at least two x and y values for interpolation.");
  }
  if (x.size() != y.size()) {
    throw std::range_error("CubicSpline requires x and y to have the same size.");
  }
  // Check that x is strictly ascending
  for (int64_t i = 1; i < x.size(); ++i) {
    if (x[i] <= x[i-1]) {
      throw std::runtime_error("CubicSpline requires x values to be strictly ascending.");
    }
  }
  
  // Generate the interpolation data
  generate_spline();
}

float CubicSpline::operator()(const float& v) const { return interpolate_spline(v); }
void CubicSpline::set_below_bounds(const float& v) { low = v; }
void CubicSpline::set_above_bounds(const float& v) { high = v; }
void CubicSpline::get_bounds(const float& value, uint64_t& lower, uint64_t& upper) const {
  upper = std::distance(points.begin(), std::lower_bound(points.begin(), points.end(), value));
  lower = upper - 1;
}

void CubicSpline::generate_spline() {
  uint64_t N = points.size();
  deriv.resize(N, 0.f);
  
  
  // Value of derivative can be found by solving tridiagonal linear equation system Ak = b (solve for k)
  std::vector<float> upper(N - 1, 0), middle(N, 0), lower(N - 1, 0), b(N, 0);
  
  // Setup the system
  upper[0] = points[1] - points[0];
  middle[0] = 2 * upper[0];
  b[0] = 3 * (values[1] - values[0]);
  for(uint64_t i = 1; i < N-1; ++i) {
    lower[i - 1] = points[i] - points[i - 1];
    upper[i] = points[i + 1] - points[i];
    middle[i] = 2 * (lower[i - 1] + upper[i]);
    b[i] = 3 * (values[i + 1] - values[i - 1]);
  }
  lower[N - 2] = points[N - 1] - points[N - 2];
  middle[N - 1] = 2 * lower[N - 2];
  b[N - 1] = 3 * (values[N - 1] - values[N - 2]);
  
  // Solve the system
  upper[0] /= middle[0];
  b[0] /= middle[0];
  for (uint64_t i = 1; i < N; ++i) {
    float tmp = 1 / std::fma(-upper[i - 1], lower[i - 1], middle[i]);
    if (i < N - 1) { upper[i] *= tmp; }
    b[i] -= b[i - 1] * lower[i - 1];
    b[i] *= tmp;
  }
  
  // Set the derivatives
  deriv[N - 1] = b[N - 1];
  for(uint64_t i = N - 2; i != std::numeric_limits<uint64_t>::max(); --i) {
    deriv[i] = std::fma(-upper[i], deriv[i + 1], b[i]);
  }
}

float CubicSpline::interpolate_spline(const float& v) const {
  if (v < points.front()) { return low; }
  if (v > points.back()) { return high; }
  
  uint64_t l, h;
  get_bounds(v, l, h);

  // q(x_val) = (1 - t(x))y1 + t(x)y2 + t(x)(1 - t(x))((1 - t(x))a + t(x)b)
  // t(x) = (x - x1) / (x2 - x1)
  // a = k1(x2 - x1) - (y2 - y1)
  // b = -k2(x2 - x1) + (y2 - y1)
  float fx = v - points[l];
  float dx = points[h] - points[l];
  float dy = values[h] - values[l];
  float t = fx / dx;
  float a = std::fma(deriv[l], dx, -dy);
  float b = std::fma(-deriv[h], dx, dy);
  float tm1 = 1 - t;
  return std::fma(tm1, values[l], std::fma(values[h], t, t * tm1 * (tm1 * a + t * b)));
  
}

}
