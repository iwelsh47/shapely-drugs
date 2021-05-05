#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <numeric>

#include <Shapely/coordinates.hpp>

namespace shapely {

using namespace Eigen;
using namespace std;

/* Constructors */
Cartesian::Cartesian() : positions(10, 3), num_points(0) { }
Cartesian::Cartesian(uint32_t count) : positions((count < 10) ? 10 : count, 3), num_points(0) { }
//Cartesian::Cartesian(const Cartesian& other) : positions(other.positions), num_points(other.num_points) { }
//Cartesian::Cartesian(Cartesian&& other) : positions(std::move(other.positions)), num_points(other.num_points) { }

/* Modifications */
Cartesian::Point Cartesian::NewPoint(float x, float y, float z) {
  CheckResize();
  positions.row(num_points) << x, y, z;
  return positions.row(num_points++);
}
void Cartesian::RemovePoint(uint32_t idx) {
  // Do nothing if idx is invalid, ie too large.
  if (idx >= num_points) { return; }
  --num_points;
  // Shift all rows below removal up one by using blocks
  positions.block(idx, 0, num_points - idx, 3) = positions.block(idx + 1, 0, num_points - idx, 3).eval();
}
void Cartesian::CheckResize(uint32_t newsize)  {
  if (!newsize && num_points == positions.rows()) {
    positions.conservativeResize(2 * positions.rows(), Eigen::NoChange_t());
  } else if (newsize && newsize > positions.rows()) {
    positions.conservativeResize(newsize, Eigen::NoChange_t());
  }
}

/* Access */
Cartesian::Point Cartesian::operator[](uint32_t idx) {
  if (idx >= num_points) { throw std::out_of_range("Requested index out of range"); }
  return positions.row(idx);
}
Cartesian::ConstPoint Cartesian::operator[](uint32_t idx) const {
  if (idx >= num_points) { throw std::out_of_range("Requested index out of range"); }
  return positions.row(idx);
}

/* Geometric operations */
Cartesian& TranslateInplace(Cartesian& pts, const Vector3f& delta, float scale) {
  Vector3f d = delta * scale;
  for (uint32_t i = 0; i < pts.NumPoints(); ++i) { pts[i] -= d; }
  return pts;
}
Cartesian Translate(const Cartesian& pts, const Vector3f& delta, float scale) {
  Cartesian result(pts);
  TranslateInplace(result, delta, scale);
  return result;
}
Cartesian& CentreInplace(Cartesian& pts, const float* weights) {
  // Calculate the current centre.
  Vector3f c;
  VectorXf w = VectorXf::Ones(pts.NumPoints());
  if (weights != nullptr) {
    memcpy(w.data(), weights, pts.NumPoints() * sizeof(float));
  }
  c.x() = pts.x().dot(w);
  c.y() = pts.y().dot(w);
  c.z() = pts.z().dot(w);
  c /= w.sum();
  
  return TranslateInplace(pts, c);
}
Cartesian Centre(const Cartesian& pts, const float* weights) {
  Cartesian result(pts);
  CentreInplace(result, weights);
  return result;
}
Cartesian& RotateInplace(Cartesian& pts, const Matrix3f& rotmat) {
  for (uint32_t i = 0; i < pts.NumPoints(); ++i) {
    pts[i] = (pts[i] * rotmat).eval();
  }
  return pts;
}
Cartesian Rotate(const Cartesian& pts, const Matrix3f& rotmat) {
  Cartesian result(pts);
  RotateInplace(result, rotmat);
  return result;
}

Matrix3f PrincipalAxesAlignment(const Cartesian& pts, const float* weights) {
  // If no points provided, return the identity matrix
  if (!pts.NumPoints()) { return Matrix3f::Identity(); }
  
  // Setup a vector for the weights
  VectorXf w = VectorXf::Ones(pts.NumPoints());
  if (weights != nullptr) {
    memcpy(w.data(), weights, pts.NumPoints() * sizeof(float));
  }
  
  // Scale the x y z columns by weights
  VectorXf x = pts.x().cwiseProduct(w);
  VectorXf y = pts.y().cwiseProduct(w);
  VectorXf z = pts.z().cwiseProduct(w);
  
  // Calculate the inertia tensor
  Matrix3f tensor;
  /* Ixx */ tensor(0,0) = (y.cwiseProduct(y) + z.cwiseProduct(z)).sum();
  /* Ixy */ tensor(0,1) = -x.cwiseProduct(y).sum();
  /* Iyy */ tensor(1,1) = (x.cwiseProduct(x) + z.cwiseProduct(z)).sum();
  /* Ixz */ tensor(0,2) = -x.cwiseProduct(z).sum();
  /* Iyz */ tensor(1,2) = -y.cwiseProduct(z).sum();
  /* Izz */ tensor(2,2) = (x.cwiseProduct(x) + y.cwiseProduct(y)).sum();
  /* Iyx */ tensor(1,0) = tensor(0,1);
  /* Izx */ tensor(2,0) = tensor(0,2);
  /* Izy */ tensor(2,1) = tensor(1,2);
  
  // Diagonalise the inertia tensor to obtain the rotation matrix
  EigenSolver<Matrix3f> eigensolver(tensor);
  Matrix3f eigensolution;
  EigenSolver<Matrix3f>::EigenvectorsType result = eigensolver.eigenvectors();
  for (uint32_t i : {0, 1, 2}) {
    for (uint32_t j : {0, 1, 2}) {
      eigensolution(i, j) = result(i, j).real();
    }
  }
  return eigensolution;
}

template <typename NumT>
inline bool approx_equal(NumT a, NumT b, NumT epsilon = std::numeric_limits<NumT>::epsilon() * 100,
                         NumT scale = 1) {
  return (std::fabs(a - b) < epsilon * (scale + std::max(std::fabs(a), std::fabs(b))));
}

Matrix3f RotationMatrix(const Eigen::Vector3f& axis, float angle) {
  // Calculate sin/cos of angle
  float sin_ = std::sin(angle);
  float cos_ = std::cos(angle);
  float cos_1 = 1 - cos_;
  
  // Extract xyz from the axis
  float x = axis.x(), y = axis.y(), z = axis.z();
  
  // Calculate the squared norm
  float sq_norm = axis.squaredNorm();
  
  // Normalise xyz if needed
  if (!approx_equal(sq_norm, 1.0f)) {
    float norm = 1 / std::sqrt(sq_norm);
    x *= norm; y *= norm; z *= norm;
  }
  Matrix3f res;
  res << cos_ + x * x * cos_1, y * x * cos_1 + z * sin_, z * x * cos_1 - y * sin_,
         x * y * cos_1 - z * sin_, cos_ + y * y * cos_1, z * y * cos_1 + x * sin_,
         x * z * cos_1 + y * sin_, y * z * cos_1 - x * sin_, cos_ + z * z * cos_1;
  return res;
}

CartesianLimits Limits(const Cartesian& pts, const float* weights) {
  CartesianLimits limit{numeric_limits<float>::max(),
                        numeric_limits<float>::max(),
                        numeric_limits<float>::max(),
                        numeric_limits<float>::lowest(),
                        numeric_limits<float>::lowest(),
                        numeric_limits<float>::lowest()};
  
  // Generate weights if weights is nullptr
  bool generate = weights == nullptr;
  if (generate) {
    weights = new float[pts.NumPoints()];
    for (size_t i = 0; i < pts.NumPoints(); ++i) { ((float*)weights)[i] = 1.0; }
  }
  
  // Calculate the min/max
  for (uint32_t i = 0; i < pts.NumPoints(); ++i) {
    limit.min_x = min(limit.min_x, pts[i].x() - weights[i]);
    limit.min_y = min(limit.min_y, pts[i].y() - weights[i]);
    limit.min_z = min(limit.min_z, pts[i].z() - weights[i]);
    
    limit.max_x = max(limit.max_x, pts[i].x() + weights[i]);
    limit.max_y = max(limit.max_y, pts[i].y() + weights[i]);
    limit.max_z = max(limit.max_z, pts[i].z() + weights[i]);
  }
  
  // Cleanup if generated weights
  if (generate) { delete [] weights; }
  
  return limit;
}

}
