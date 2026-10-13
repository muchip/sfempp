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

#ifndef GEOMETRY_LSHAPE_H_
#define GEOMETRY_LSHAPE_H_

namespace SFEMPP {
struct LshapeCurve {
  LshapeCurve(unsigned int P) : _P(P) {}
  Eigen::Vector2d f(double s) const {
    Eigen::Vector2d retVal;
    switch (_P) {
      case 0:
        retVal << 0.5 + s, -0.5;
        return retVal;
      case 1:
        retVal << 1.5, -0.5 + s;
        return retVal;
      case 2:
        retVal << 1.5 - 2 * s, 0.5;
        return retVal;
      case 3:
        retVal << -0.5, 0.5 - 2 * s;
        return retVal;
      case 4:
        retVal << -0.5 + s, -1.5;
        return retVal;
      case 5:
        retVal << 0.5, -1.5 + s;
        return retVal;

      default:
        return 0 * retVal;
    }
  }
  Eigen::Vector2d df(double s) const {
    Eigen::Vector2d retVal;
    switch (_P) {
      case 0:
        retVal << 1, 0;
        return retVal;
      case 1:
        retVal << 0, 1;
        return retVal;
      case 2:
        retVal << -2, 0;
        return retVal;
      case 3:
        retVal << 0, -2;
        return retVal;
      case 4:
        retVal << 1, 0;
        return retVal;
      case 5:
        retVal << 0, 1;
        return retVal;

      default:
        return 0 * retVal;
    }
  }
  static constexpr unsigned int numPatches = 6;
  unsigned int _P;
};

typedef Parametrisation<LshapeCurve::numPatches, LshapeCurve> Lshape;
}  // namespace SFEMPP
#endif
