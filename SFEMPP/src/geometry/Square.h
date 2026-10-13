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

#ifndef GEOMETRY_SQUARE_H_
#define GEOMETRY_SQUARE_H_

namespace SFEMPP {
struct SquareCurve {
  SquareCurve(unsigned int P) : _P(P) {}
  Eigen::Vector2d f(double s) const {
    Eigen::Vector2d retVal;
    switch (_P) {
      case 0:
        retVal << 0.5, -0.5 + s;
        return retVal;
      case 1:
        retVal << 0.5 - s, 0.5;
        return retVal;
      case 2:
        retVal << -0.5, 0.5 - s;
        return retVal;
      case 3:
        retVal << -0.5 + s, -0.5;
        return retVal;
      default:
        return 0 * retVal;
    }
  }
  Eigen::Vector2d df(double s) const {
    Eigen::Vector2d retVal;
    switch (_P) {
      case 0:
        retVal << 0, 1;
        return retVal;
      case 1:
        retVal << -1, 0;
        return retVal;
      case 2:
        retVal << 0, -1;
        return retVal;
      case 3:
        retVal << 1, 0;
        return retVal;
      default:
        return 0 * retVal;
    }
  }
  static constexpr unsigned int numPatches = 4;
  unsigned int _P;
};

typedef Parametrisation<SquareCurve::numPatches, SquareCurve> Square;
}  // namespace SFEMPP
#endif
