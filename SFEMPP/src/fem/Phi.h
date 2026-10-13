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

#ifndef SFEMPP_FEM_PHI_H_
#define SFEMPP_FEM_PHI_H_

namespace SFEMPP {
namespace FEM {
/*
 *  class that implements the shape functions on the reference element
 */
class Phi {
 public:
  double f(int i, const Eigen::Vector2d &a) const {
    switch (i) {
      case 0:
        return phi0(a);
      case 1:
        return phi1(a);
      case 2:
        return phi2(a);
      default:
        std::cerr << "there are only 3 basis functions (0,1,2).\n";
        return 0. / 0.;
    }
  }
  double dx(int i, const Eigen::Vector2d &a) const {
    switch (i) {
      case 0:
        return dphi0dx(a);
      case 1:
        return dphi1dx(a);
      case 2:
        return dphi2dx(a);
      default:
        std::cerr << "there are only 3 basis functions (0,1,2).\n";
        return 0. / 0.;
    }
  }
  double dy(int i, const Eigen::Vector2d &a) const {
    switch (i) {
      case 0:
        return dphi0dy(a);
      case 1:
        return dphi1dy(a);
      case 2:
        return dphi2dy(a);
      default:
        std::cerr << "there are only 3 basis functions (0,1,2).\n";
        return 0. / 0.;
    }
  }

 private:
  double phi0(const Eigen::Vector2d &a) const { return 1 - a(0) - a(1); }
  double phi1(const Eigen::Vector2d &a) const { return a(0); }
  double phi2(const Eigen::Vector2d &a) const { return a(1); }
  double dphi0dx(const Eigen::Vector2d &a) const {
    (void)a;
    return -1;
  }
  double dphi1dx(const Eigen::Vector2d &a) const {
    (void)a;
    return 1;
  }
  double dphi2dx(const Eigen::Vector2d &a) const {
    (void)a;
    return 0;
  }
  double dphi0dy(const Eigen::Vector2d &a) const {
    (void)a;
    return -1;
  }
  double dphi1dy(const Eigen::Vector2d &a) const {
    (void)a;
    return 0;
  }
  double dphi2dy(const Eigen::Vector2d &a) const {
    (void)a;
    return 1;
  }
};
}  // namespace FEM
}  // namespace SFEMPP
#endif
