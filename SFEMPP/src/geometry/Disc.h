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

#ifndef GEOMETRY_DISC_H_
#define GEOMETRY_DISC_H_

namespace SFEMPP {
struct Circle {
  Circle(unsigned int P) : _P(P) {}
  Eigen::Vector2d f(double s) const {
    double arg = s + double(_P);
    double c = 2 * M_PI / double(numPatches);
    Eigen::Vector2d retVal;
    retVal << std::cos(c * arg), std::sin(c * arg);
    return retVal;
  }
  Eigen::Vector2d df(double s) const {
    double arg = s + double(_P);
    double c = 2 * M_PI / double(numPatches);
    Eigen::Vector2d retVal;
    retVal << -c * std::sin(c * arg), c * std::cos(c * arg);
    return retVal;
  }
  static constexpr unsigned int numPatches = 4;
  unsigned int _P;
};

typedef Parametrisation<Circle::numPatches, Circle> Disc;
}  // namespace SFEMPP
#endif
