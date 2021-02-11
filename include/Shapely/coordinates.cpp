#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <Shapely/coordinates.hpp>

namespace shapely {

using namespace Eigen;

/* Constructors */
Cartesian::Cartesian() : positions(10, 3), num_points(0) { }
Cartesian::Cartesian(uint32_t count) : positions((count < 10) ? 10 : count, 3), num_points(0) { }
Cartesian::Cartesian(const Cartesian& other) : positions(other.positions), num_points(other.num_points) { }
Cartesian::Cartesian(Cartesian&& other) : positions(std::move(other.positions)), num_points(other.num_points) { }

/* Modifications */
Cartesian::Point Cartesian::NewPoint(double x, double y, double z) {
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
Cartesian& TranslateInplace(Cartesian& pts, const Vector3d& delta, double scale) {
  Vector3d d = delta * scale;
  for (uint32_t i = 0; i < pts.NumPoints(); ++i) { pts[i] -= d; }
  return pts;
}
Cartesian Translate(const Cartesian& pts, const Vector3d& delta, double scale) {
  Cartesian result(pts);
  TranslateInplace(result, delta, scale);
  return result;
}
Cartesian& CentreInplace(Cartesian& pts, const double* weights) {
  // Calculate the current centre.
  Vector3d c;
  VectorXd w = VectorXd::Ones(pts.NumPoints());
  if (weights != nullptr) {
    memcpy(w.data(), weights, pts.NumPoints() * sizeof(double));
  }
  c.x() = pts.x().dot(w);
  c.y() = pts.y().dot(w);
  c.z() = pts.z().dot(w);
  c /= w.sum();
  
  return TranslateInplace(pts, c);
}
Cartesian Centre(const Cartesian& pts, const double* weights) {
  Cartesian result(pts);
  CentreInplace(result, weights);
  return result;
}
Cartesian& RotateInplace(Cartesian& pts, const Matrix3d& rotmat) {
  for (uint32_t i = 0; i < pts.NumPoints(); ++i) {
    pts[i] = (pts[i] * rotmat).eval();
  }
  return pts;
}
Cartesian Rotate(const Cartesian& pts, const Matrix3d& rotmat) {
  Cartesian result(pts);
  RotateInplace(result, rotmat);
  return result;
}

Matrix3d PrincipalAxesAlignment(const Cartesian& pts, const double* weights) {
  // If no points provided, return the identity matrix
  if (!pts.NumPoints()) { return Matrix3d::Identity(); }
  
  // Setup a vector for the weights
  VectorXd w = VectorXd::Ones(pts.NumPoints());
  if (weights != nullptr) {
    memcpy(w.data(), weights, pts.NumPoints() * sizeof(double));
  }
  
  // Scale the x y z columns by weights
  VectorXd x = pts.x().cwiseProduct(w);
  VectorXd y = pts.y().cwiseProduct(w);
  VectorXd z = pts.z().cwiseProduct(w);
  
  // Calculate the inertia tensor
  Matrix3d tensor;
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
  EigenSolver<Matrix3d> eigensolver(tensor);
  Matrix3d eigensolution;
  EigenSolver<Matrix3d>::EigenvectorsType result = eigensolver.eigenvectors();
  for (uint32_t i : {0, 1, 2}) {
    for (uint32_t j : {0, 1, 2}) {
      eigensolution(i, j) = result(i, j).real();
    }
  }
  return eigensolution;
}

}
