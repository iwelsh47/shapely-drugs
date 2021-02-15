#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <cstdint>
#include <iostream>

namespace shapely {

/** \brief Container used to store cartesian coordinates.
 *  \details Internally, uses an Eigen::Matrix to store and manipulate the coordinates
 */
class Cartesian {
  using Matrix = Eigen::MatrixX3d;
public:
  /** Type for interacting with the underlying storage. */
  using Coordinates = typename Matrix::FixedBlockXpr<Eigen::Dynamic, 3>::Type;
  using ConstCoordinates = typename Matrix::ConstFixedBlockXpr<Eigen::Dynamic, 3>::Type;
  
  /** Type for interacting with a single column (eg x values) of the storage. */
  using Axis = typename Matrix::FixedBlockXpr<Eigen::Dynamic, 1>::Type;
  using ConstAxis = typename Matrix::ConstFixedBlockXpr<Eigen::Dynamic, 1>::Type;
  
  /** Type for interacting with a single row of the storage. */
  using Point = typename Matrix::FixedBlockXpr<1, 3>::Type;
  using ConstPoint = typename Matrix::ConstFixedBlockXpr<1, 3>::Type;
  
  /** \brief Default constructor.
   *  \details Constructs an empty container, with space for 10 points reserved.
   */
  Cartesian();
  
  /** \brief Construct the container with space for \p count points.
   *  \details Reverses space for at least \p count points, though no points are initalised.
   *  \param count number of points to reserve space for.
   */
  Cartesian(uint32_t count);
  
  /** \brief Copy constructor.
   *  \details Constructs as a copy of another existing Cartesian instance.
   *  \param other the existing instance to copy from.
   */
  Cartesian(const Cartesian& other);
  
  /** \brief Move constructor.
   *  \details Constructs using move semantics.
   *  \param other the existing instance to construct from
   */
  Cartesian(Cartesian&& other);
  
  /** Destructor */
  ~Cartesian() { }
  
  /** \brief Add a new point to the set.
   *  \details Creates a new point with the given coordinate values. If the storage matrix does not have enough space, the reserved space is doubled in size.
   *  \param x,y,z cartesian coordinates of the point to add.
   */
  Point NewPoint(double x, double y, double z);
  inline Point NewPoint(const Eigen::RowVector3d& pt) { return NewPoint(pt.x(), pt.y(), pt.z()); }
  inline Point NewPoint(const Eigen::Vector3d& pt) { return NewPoint(pt.x(), pt.y(), pt.z()); }
  
  /** \brief Remove a point from the set.
   *  \details Removes the point at the given index from the set.
   *  \param idx point index to remove.
   */
  void RemovePoint(uint32_t idx);
  
  /** Reserve space for up to \p num points. */
  inline void Reserve(uint32_t num) { CheckResize(num); }
  
  /** Access the point at the given index. */
  Point operator[](uint32_t idx);
  ConstPoint operator[](uint32_t idx) const;
  
  /** Get the number of points in the container. */
  inline uint32_t NumPoints() const { return num_points; }
  /** Get the underlying storage matrix. */
  inline Coordinates GetCoordinates() { return positions.block<Eigen::Dynamic, 3>(0, 0, num_points, 3); }
  inline ConstCoordinates GetCoordinates() const { return positions.block<Eigen::Dynamic, 3>(0, 0, num_points, 3); }
  
  /** Access the indivdual axis set of coordinates. */
  inline Axis x() { return positions.block<Eigen::Dynamic, 1>(0, 0, num_points, 1); }
  inline Axis y() { return positions.block<Eigen::Dynamic, 1>(0, 1, num_points, 1); }
  inline Axis z() { return positions.block<Eigen::Dynamic, 1>(0, 2, num_points, 1); }
  inline ConstAxis x() const { return positions.block<Eigen::Dynamic, 1>(0, 0, num_points, 1); }
  inline ConstAxis y() const { return positions.block<Eigen::Dynamic, 1>(0, 1, num_points, 1); }
  inline ConstAxis z() const { return positions.block<Eigen::Dynamic, 1>(0, 2, num_points, 1); }
  
private:
  /** Check if resizing is required and do it if needed. */
  void CheckResize(uint32_t newsize = 0);
  
public:
  Matrix positions;
  uint32_t num_points;
};

/** Print out matrix to an ostream. */
static inline std::ostream& operator<<(std::ostream& os, const Cartesian& coords) {
  os << coords.GetCoordinates();
  return os;
}

/** Translate the points in \p pts by the given vector. Translation is performed by subtraction. */
Cartesian& TranslateInplace(Cartesian& pts, const Eigen::Vector3d& delta, double scale = 1.0);
Cartesian Translate(const Cartesian& pts, const Eigen::Vector3d& delta, double scale = 1.0);
inline Cartesian& TranslateInplace(Cartesian& pts, double dx, double dy, double dz, double scale = 1.0) {
  return TranslateInplace(pts, Eigen::Vector3d(dx, dy, dz), scale);
}
inline Cartesian Translate(const Cartesian& pts, double dx, double dy, double dz, double scale = 1.0) {
  return Translate(pts, Eigen::Vector3d(dx, dy, dz), scale);
}

/** Centre the coordinates at (0, 0, 0), optionally with weights. */
Cartesian& CentreInplace(Cartesian& pts, const double* weights = nullptr);
Cartesian Centre(const Cartesian& pts, const double* weights = nullptr);

/** Rotate the coordinates with the given rotation matrix. */
Cartesian& RotateInplace(Cartesian& pts, const Eigen::Matrix3d& rotmat);
Cartesian Rotate(const Cartesian& pts, const Eigen::Matrix3d& rotmat);

/** Calculate the rotation matrix to obtain principal axes alignment. */
Eigen::Matrix3d PrincipalAxesAlignment(const Cartesian& pts, const double* weights = nullptr);
/** Calculate the rotation matrix given an axis and rotation angle. Angle in radians*/
Eigen::Matrix3d RotationMatrix(const Eigen::Vector3d& axis, double angle);

}
