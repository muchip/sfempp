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

#ifndef GEOMETRY_STAR_H_
#define GEOMETRY_STAR_H_

namespace SFEMPP {
struct StarCurve {
  StarCurve(unsigned int P) : _P(P) {}
  Eigen::Vector2d f(double s) const {
    double arg = s + double(_P);
    double c = 2 * M_PI / double(numPatches);
    double rad = 1 + 0.25 * std::sin(5 * c * arg);
    Eigen::Vector2d retVal;
    retVal << rad * std::cos(c * arg), rad * std::sin(c * arg);
    return retVal;
  }
  Eigen::Vector2d df(double s) const {
    double arg = s + double(_P);
    double c = 2 * M_PI / double(numPatches);
    double rad = 1 + 0.25 * std::sin(5 * c * arg);
    double drad = 1.25 * c * std::cos(5 * c * arg);
    Eigen::Vector2d retVal;
    retVal << -rad * c * std::sin(c * arg) + drad * std::cos(c * arg),
        rad * c * std::cos(c * arg) + drad * std::sin(c * arg);
    return retVal;
  }
  static constexpr unsigned int numPatches = 7;
  unsigned int _P;
};

typedef Parametrisation<StarCurve::numPatches, StarCurve> Star;
}  // namespace SFEMPP
#endif
