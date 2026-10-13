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

#ifndef SFEMPP_FEM_INTERPOLATION_H_
#define SFEMPP_FEM_INTERPOLATION_H_

namespace SFEMPP {
namespace MMG {
template <typename T>
Eigen::VectorXd interpolation(const Mesh& msh, const T& func) {
  auto PC = msh.PC();
  Eigen::VectorXd retVal(PC.cols());
  for (auto i = 0; i < PC.cols(); ++i) retVal(i) = func(PC.col(i));
  return retVal;
}

Eigen::VectorXd projection(const Mesh& msh, const Eigen::VectorXd& f, int m1,
                           int m2) {
  Eigen::VectorXd retVal = f;
  auto maxLvl = msh.maxlvl();
  eigen_assert(m2 <= maxLvl);
  for (auto i = m1 + 1; i <= m2; ++i) prolongation(retVal, msh, i);
  return retVal;
}
}  // namespace MMG
}  // namespace SFEMPP
#endif
