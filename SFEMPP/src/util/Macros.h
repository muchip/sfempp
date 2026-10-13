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

#ifndef SFEMPP_UTIL_MACROS_H_
#define SFEMPP_UTIL_MACROS_H_

#include <limits>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841971693993751
#endif

namespace SFEMPP {
// define primitive types used throughout the toolbox
#define SFEMPP_INDEX unsigned int
#define SFEMPP_SCALAR double
#define SFEMPP_INF std::numeric_limits<SFEMPP_SCALAR>::infinity()
#define SFEMPP_ZERO_TOLERANCE std::numeric_limits<SFEMPP_SCALAR>::epsilon()

typedef SFEMPP_INDEX Index;

typedef SFEMPP_SCALAR Scalar;

}  // namespace SFEMPP

#endif
