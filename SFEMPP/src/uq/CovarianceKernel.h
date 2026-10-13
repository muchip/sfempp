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
#ifndef SFEMPP_UQ_COVKERNEL_H_
#define SFEMPP_UQ_COVKERNEL_H_

namespace SFEMPP {
class CovarianceKernelMF {
 public:
  CovarianceKernelMF(const Eigen::MatrixXd &pts)
      : pts_(pts), dim_(2 * pts.cols()){};
  unsigned int dim() const { return dim_; }
  double operator()(int i, int j) const {
    // const double l = 0.5;
    int np = pts_.cols();
    // case indexes matrix according to [A_0 A_1 ; A_2 A_3]
    int CASE = 2 * (i / np) + (j / np);
    // now crop i and j to sizes smaller than np
    i = i % np;
    j = j % np;
    const double r = (pts_.col(i) - pts_.col(j)).norm();
    const double normx = pts_.col(i).norm();
    const double normy = pts_.col(j).norm();
    const double kernel = 9./16. * exp(-0.01 * r * r);
    const Eigen::Vector2d dirx = pts_.col(i);
    const Eigen::Vector2d diry = pts_.col(j);
    // incorporate variance 1 here by adding a factor 3 to the kernel
    switch (CASE) {
      case 0:
        return kernel;
        break;
      case 1:
        return 0.;
        break;
      case 2:
        return 0.;
        break;
      case 3:
        return 0.;
        break;
    }
    return SFEMPP_INF;
  }

 private:
  const Eigen::MatrixXd &pts_;
  unsigned int dim_;
};
}  // namespace SFEMPP
#endif
