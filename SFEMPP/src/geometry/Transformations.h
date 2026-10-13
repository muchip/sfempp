// This file is part of SFEMPP, a lightwight stochastic FEM package.
//
// Copyright (c) 2023, Michael Multerer
//
// All rights reserved.
//
// This source code is subject to the GNU Affero General Public License v3.0
// license and without any warranty, see <https://github.com/muchip/FMCA>
// for further information.
//

#ifndef GEOMETRY_TRANSFORMATIONS_H_
#define GEOMETRY_TRANSFORMATIONS_H_

namespace SFEMPP {
namespace Transformations {
struct Curve {
  Curve(Eigen::Vector2d (*f)(double), Eigen::Vector2d (*df)(double)) {
    this->f = f;
    this->df = df;
  }
  Eigen::Vector2d (*f)(double);
  Eigen::Vector2d (*df)(double);
};

/**
 *  functions related to the affine mapping from the reference
 *  triangle to the geometry
 */
// affine mapping
static Eigen::Vector2d F_T(const Eigen::Vector2d &a, const Eigen::Vector2d &p0,
                           const Eigen::Vector2d &p1,
                           const Eigen::Vector2d &p2) {
  return p0 + a(0) * (p1 - p0) + a(1) * (p2 - p0);
}
// derivative in x direction
static Eigen::Vector2d dF_Tdx(const Eigen::Vector2d &a,
                              const Eigen::Vector2d &p0,
                              const Eigen::Vector2d &p1,
                              const Eigen::Vector2d &p2) {
  (void)a;
  (void)p2;
  return p1 - p0;
}
// derivative in y direction
static Eigen::Vector2d dF_Tdy(const Eigen::Vector2d &a,
                              const Eigen::Vector2d &p0,
                              const Eigen::Vector2d &p1,
                              const Eigen::Vector2d &p2) {
  (void)a;
  (void)p1;
  return p2 - p0;
}
// functional determinant of the element
static double det_dF_T(Eigen::Vector2d a, Eigen::Vector2d p0,
                       Eigen::Vector2d p1, Eigen::Vector2d p2) {
  (void)a;
  return ((p1(0) - p0(0)) * (p2(1) - p0(1)) -
          (p1(1) - p0(1)) * (p2(0) - p0(0)));
}
// inverse mapping from the triangle p0,p1,p2 to the reference element
static Eigen::Vector2d invF_T(const Eigen::Vector2d &a,
                              const Eigen::Vector2d &p0,
                              const Eigen::Vector2d &p1,
                              const Eigen::Vector2d &p2) {
  Eigen::Matrix2d invA;
  double invDet = 1. / det_dF_T(a, p0, p1, p2);
  invA << (p2(1) - p0(1)) * invDet, -(p2(0) - p0(0)) * invDet,
      -(p1(1) - p0(1)) * invDet, (p1(0) - p0(0)) * invDet;
  return invA * (a - p0);
}
/**
 *  functions related to the affine mapping from the reference
 *  triangle to the geometry using Zenizek triangles
 *  here, the connecting edge (p1,p2) is given by a curve gamma from [0,1]
 * with gamma.f(0) = p1 and gamma.f(1) = p2
 */
// mapping from the reference element to the triangle p0,p1,p2
template <typename T>
static Eigen::Vector2d tildeF_T(const Eigen::Vector2d &a,
                                const Eigen::Vector2d &p0, const T &gamma) {
  if (std::abs(a(1) - 1) < 2 * SFEMPP_ZERO_TOLERANCE) return gamma.f(1);
  if (std::abs(a(0) - 1) < 2 * SFEMPP_ZERO_TOLERANCE) return gamma.f(0);
  if (std::abs(a(0)) + std::abs(a(1)) < 2 * SFEMPP_ZERO_TOLERANCE)
    return p0;
  else
    return (1 - a(0) - a(1)) * (p0 + a(1) / (1 - a(1)) * gamma.f(1)) +
           a(0) / (1 - a(1)) * gamma.f(a(1));
}
// derivative in x direction
template <typename T>
static Eigen::Vector2d dtildeF_Tdx(const Eigen::Vector2d &a,
                                   const Eigen::Vector2d &p0, const T &gamma) {
  if (std::abs(a(1) - 1) < 2 * SFEMPP_ZERO_TOLERANCE)
    return gamma.f(1) - p0 - gamma.df(1);
  else
    return 1 / (1 - a(1)) * gamma.f(a(1)) -
           (p0 + a(1) / (1 - a(1)) * gamma.f(1));
}

// derivative in y direction
template <typename T>
static Eigen::Vector2d dtildeF_Tdy(const Eigen::Vector2d &a,
                                   const Eigen::Vector2d &p0, const T &gamma) {
  if (std::abs(a(1) - 1) < 2 * SFEMPP_ZERO_TOLERANCE)
    return gamma.f(1) - p0;
  else
    return (1 - a(0) / ((1 - a(1)) * (1 - a(1)))) * gamma.f(1) - p0 +
           a(0) / ((1 - a(1)) * (1 - a(1))) * gamma.f(a(1)) +
           a(0) / (1 - a(1)) * gamma.df(a(1));
}
// functional determinant of the element
template <typename T>
static double det_dtildeF_T(const Eigen::Vector2d &a, const Eigen::Vector2d &p0,
                            const T &gamma) {
  Eigen::Vector2d v1 = dtildeF_Tdx(a, p0, gamma);
  Eigen::Vector2d v2 = dtildeF_Tdy(a, p0, gamma);
  return v1(0) * v2(1) - v1(1) * v2(0);
}
}  // namespace Transformations
}  // namespace SFEMPP
#endif
