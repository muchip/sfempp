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

#ifndef SFEMPP_FEM_RIGHTHANDSIDE_H_
#define SFEMPP_FEM_RIGHTHANDSIDE_H_

namespace SFEMPP {
namespace FEM {
template <typename T>
Eigen::VectorXd computeLoadVector(const Mesh &msh, const T &func) {
  using namespace Transformations;
  auto P = msh.PC();
  auto nP = P.cols();
  auto F = msh.Ffull().back();
  auto B = msh.B();
  Phi phi;
  Quadrature Q;
  Eigen::VectorXd load(nP);
  load.setZero();
  for (auto i = 0; i < F.cols(); ++i) {
    Eigen::Vector3d ints;
    ints.setZero();
    // determine all integrals on the element in question
    for (auto j = 0; j < Q._w.size(); ++j) {
      Eigen::Vector2d xi = Q._xi.col(j);
      double w = Q._w(j) *
                 det_dF_T(xi, P.col(F(0, i)), P.col(F(1, i)), P.col(F(2, i))) *
                 func(F_T(xi, P.col(F(0, i)), P.col(F(1, i)), P.col(F(2, i))));
      for (auto k = 0; k < 3; ++k) ints(k) += w * phi.f(k, xi);
    }
    // write integrals to particular entries
    for (auto j = 0; j < 3; ++j)
      if (B(F(j, i)) != 0) load(F(j, i)) += ints(j);
  }
  return load;
}

}  // namespace FEM
}  // namespace SFEMPP
#endif
