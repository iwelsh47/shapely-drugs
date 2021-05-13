#pragma once

#include <vector>

namespace shapely {

class CubicSpline {
private:
  void generate_spline();
  float interpolate_spline(const float& value) const;
  void get_bounds(const float& value, uint64_t& lower, uint64_t& upper) const;
  
public:
  CubicSpline() = default;
  CubicSpline(const std::vector<float>& x_values, const std::vector<float>& y_values);
  ~CubicSpline() = default;
  
  float operator()(const float& v) const;
  
  void set_below_bounds(const float& new_v);
  void set_above_bounds(const float& new_v);
  
  inline const std::vector<float>& get_interpolation_points() const { return points; }
  inline const std::vector<float>& get_interpolation_values() const { return values; }
  inline const std::vector<float>& get_derivatives() const { return deriv; }
  
private:
  std::vector<float> points, values, deriv;
  float low, high;
};

}
