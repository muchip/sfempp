// This file is part of SFEMPP, a lightwight stochastic FEM package.
//
// Copyright (c) 2023, Michael Multerer
//
// All rights reserved.
//
// This source code is subject to the GNU Affero General Public License v3.0
// license and without any warranty, see <https://github.com/muchip/SFEMPP>
// for further information.
//

#ifndef SFEMPP_FEM_NORMS_H_
#define SFEMPP_FEM_NORMS_H_

namespace SFEMPP {
namespace FEM {

/*
 *  determines the difference of the gradient of f and the gradient
 *  given by func
 */
template <typename __T>
double computeH1error(const Mesh &msh, const Eigen::VectorXd &f,
                      const __T &func) {
  using namespace Transformations;
  auto P = msh.PC();
  auto nP = msh.nP();
  auto F = msh.Ffull();
  auto B = msh.B();
  Eigen::Vector2d gradf;
  Eigen::Matrix2d invJ;
  Eigen::Vector2d gradPhi;
  Phi phi;
  Quadrature Q;
  double normf = 0;
  for (auto i = 0; i < F.back().cols(); ++i) {
    Eigen::Vector3i p = F.back().col(i);
    for (auto j = 0; j < Q._w.size(); ++j) {
      Eigen::Vector2d xi = Q._xi.col(j);
      Eigen::Vector2d FTxi = F_T(xi, P.col(p(0)), P.col(p(1)), P.col(p(2)));
      gradf.setZero();
      auto a = dF_Tdx(xi, P.col(p(0)), P.col(p(1)), P.col(p(2)));
      auto b = dF_Tdy(xi, P.col(p(0)), P.col(p(1)), P.col(p(2)));
      double c = det_dF_T(xi, P.col(p(0)), P.col(p(1)), P.col(p(2)));
      invJ << b(1) / c, -b(0) / c, -a(1) / c, a(0) / c;
      for (auto k = 0; k < 3; ++k) {
        gradPhi << phi.dx(k, xi), phi.dy(k, xi);
        gradf += gradPhi * f(p(k));
      }
      // repair scaling of invers by dividing with detJ
      normf +=
          (invJ.transpose() * gradf - func(FTxi)).squaredNorm() * Q._w(j) * c;
    }
  }
  return std::sqrt(normf);
}
/*
 *  determines the L2 norm of a finite element function on the mesh msh
 */
template <typename __T>
double computeL2error(const Mesh &msh, const Eigen::VectorXd &f,
                      const __T &func) {
  using namespace Transformations;
  auto P = msh.PC();
  auto nP = msh.nP();
  auto F = msh.Ffull();
  auto B = msh.B();
  Phi phi;
  Quadrature Q;
  double normf = 0;
  for (auto i = 0; i < F.back().cols(); ++i) {
    Eigen::Vector3i p = F.back().col(i);
    for (auto j = 0; j < Q._w.size(); ++j) {
      Eigen::Vector2d xi = Q._xi.col(j);
      Eigen::Vector2d FTxi = F_T(xi, P.col(p(0)), P.col(p(1)), P.col(p(2)));
      double c = 0;
      for (auto k = 0; k < 3; ++k) c += phi.f(k, xi) * f(p(k));
      normf += std::pow(c - func(FTxi), 2.) *
               det_dF_T(xi, P.col(p(0)), P.col(p(1)), P.col(p(2)));
    }
  }
  return std::sqrt(normf);
}

/*
 *  determines the L2 norm of a finite element function on the mesh msh
 */
double computeL2norm(const Eigen::MatrixXd &P, const Eigen::MatrixXi &E,
                     const Eigen::VectorXd &f) {
  using namespace Transformations;
  Phi phi;
  Quadrature Q;
  double normf = 0;
  for (auto i = 0; i < E.cols(); ++i) {
    Eigen::Vector3i p = E.col(i);
    for (auto j = 0; j < Q._w.size(); ++j) {
      Eigen::Vector2d xi = Q._xi.col(j);
      double c = 0;
      for (auto k = 0; k < 3; ++k) c += phi.f(k, xi) * f(p(k));
      normf += c * c * Q._w(j) *
               abs(det_dF_T(xi, P.col(p(0)), P.col(p(1)), P.col(p(2))));
    }
  }
  return std::sqrt(normf);
}

}  // namespace FEM
}  // namespace SFEMPP
#endif
