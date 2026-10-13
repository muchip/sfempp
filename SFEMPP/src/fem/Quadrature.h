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

#ifndef SFEMPP_FEM_QUADRATURE_H_
#define SFEMPP_FEM_QUADRATURE_H_

namespace SFEMPP {
namespace FEM {
struct Quadrature {
  Quadrature(void) {
    _w.resize(7);
    _xi.resize(2, 7);
    _w << 9. / 80, (155 + sqrt(15)) / 2400, (155 + sqrt(15)) / 2400,
        (155 + sqrt(15)) / 2400, (155 - sqrt(15)) / 2400,
        (155 - sqrt(15)) / 2400, (155 - sqrt(15)) / 2400;
    _xi.row(0) << 1. / 3, (6 + sqrt(15)) / 21, (9 - 2 * sqrt(15)) / 21,
        (6 + sqrt(15)) / 21, (6 - sqrt(15)) / 21, (9 + 2 * sqrt(15)) / 21,
        (6 - sqrt(15)) / 21;
    _xi.row(1) << 1. / 3, (6 + sqrt(15)) / 21, (6 + sqrt(15)) / 21,
        (9 - 2 * sqrt(15)) / 21, (6 - sqrt(15)) / 21, (6 - sqrt(15)) / 21,
        (9 + 2 * sqrt(15)) / 21;
  }
  Eigen::MatrixXd _xi;
  Eigen::VectorXd _w;
};

struct Quadrature1D {
  Quadrature1D(void) {
    _w.resize(7);
    _xi.resize(7);
    _w << 6.47424830844348466353e-02, 1.39852695744638333951e-01,
        1.90915025252559472475e-01, 2.08979591836734693878e-01,
        1.90915025252559472475e-01, 1.39852695744638333951e-01,
        6.47424830844348466353e-02;
    _xi << 2.54460438286207377369e-02, 1.29234407200302780068e-01,
        2.97077424311301416547e-01, 5.00000000000000000000e-01,
        7.02922575688698583453e-01, 8.70765592799697219932e-01,
        9.74553956171379262263e-01;
  }
  Eigen::VectorXd _xi;
  Eigen::VectorXd _w;
};

}  // namespace FEM
}  // namespace SFEMPP
#endif
