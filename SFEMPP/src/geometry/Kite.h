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

#ifndef GEOMETRY_KITE_H_
#define GEOMETRY_KITE_H_

namespace SFEMPP {
struct KiteCurve {
  KiteCurve(unsigned int P) : _P(P) {}
  Eigen::Vector2d f(double s) const {
    double arg = s + double(_P);
    double c = 2 * M_PI / double(numPatches);
    Eigen::Vector2d retVal;
    retVal << 5 * std::cos(c * arg) - 3.25 * std::cos(2 * c * arg) + 3.25,
        7.5 * std::sin(c * arg);
    return retVal;
  }
  Eigen::Vector2d df(double s) const {
    double arg = s + double(_P);
    double c = 2 * M_PI / double(numPatches);
    Eigen::Vector2d retVal;
    retVal << -5 * c * std::sin(c * arg) + 6.5 * c * std::sin(2 * c * arg),
        7.5 * c * std::cos(c * arg);
    return retVal;
  }
  static constexpr unsigned int numPatches = 5;
  unsigned int _P;
};

typedef Parametrisation<KiteCurve::numPatches, KiteCurve> Kite;
}  // namespace SFEMPP
#endif
