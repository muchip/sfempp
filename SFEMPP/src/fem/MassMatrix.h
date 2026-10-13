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

#ifndef SFEMPP_FEM_MASSMATRIX_H_
#define SFEMPP_FEM_MASSMATRIX_H_

namespace SFEMPP {
namespace FEM {

Eigen::VectorXd intMass(const Eigen::MatrixXd &P, const Eigen::VectorXi &F,
                        const Phi &phi, const Quadrature &Q) {
  using namespace Transformations;
  Eigen::VectorXd mass(9);
  mass.setZero();
  Eigen::Vector3i p = F;
  for (auto i = 0; i < Q._w.size(); ++i) {
    Eigen::Vector2d xi = Q._xi.col(i);
    double w = Q._w(i) * det_dF_T(xi, P.col(p(0)), P.col(p(1)), P.col(p(2)));
    for (auto k = 0; k < 9; ++k)
      mass(k) += w * phi.f(k / 3, xi) * phi.f(k % 3, xi);
  }
  return mass;
}

std::vector<Eigen::SparseMatrix<double>> computeMassMatrixMG(
    const Mesh &msh, const std::vector<Eigen::MatrixXd> &massWeights,
    int l = -1) {
  auto maxLvl = msh.maxlvl();
  auto P = msh.PC();
  auto nP = msh.nP();
  auto F = msh.Ffull();
  auto B = msh.B();
  Eigen::Vector3i p;
  if (l == -1 || l > maxLvl) l = maxLvl;
  std::vector<Eigen::SparseMatrix<double>> massMatrix(l + 1);
  auto itA = massMatrix.rbegin();
  auto itF = F.rbegin();
  auto itn = nP.rbegin();
  auto itW = massWeights.rbegin();
  std::advance(itF, maxLvl - l);
  std::advance(itn, maxLvl - l);
  std::advance(itW, maxLvl - l);
  for (; itF != F.rend(); ++itA, ++itF, ++itn, ++itW) {
    Eigen::SparseMatrix<double> M(*itn, *itn);
    M.reserve((*itn) * (msh.computeMaxDegree() + 1));
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve((*itn) * (msh.computeMaxDegree() + 1) * 9);
    // add entries into sparse matrix
    for (auto i = 0; i < itF->cols(); ++i) {
      p = itF->col(i);
      for (auto j = 0; j < 3; ++j)
        for (auto k = 0; k < 3; ++k)
          if (B(p(j)) && B(p(k)))
            tripletList.push_back(
                Eigen::Triplet<double>(p(j), p(k), (*itW)(j + 3 * k, i)));
    }
    M.setFromTriplets(tripletList.begin(), tripletList.end());
    *itA = std::move(M);
  }
  return massMatrix;
}

Eigen::SparseMatrix<double> computeMassMatrix(const Mesh &msh) {
  using namespace Transformations;
  auto P = msh.PC();
  auto nP = msh.nP();
  auto F = msh.Ffull();
  auto B = msh.B();
  Phi phi;
  Quadrature Q;
  Eigen::SparseMatrix<double> M(nP.back(), nP.back());
  M.reserve((nP.back()) * (msh.computeMaxDegree() + 1));
  std::vector<Eigen::Triplet<double>> tripletList;
  tripletList.reserve((nP.back()) * (msh.computeMaxDegree() + 1) * 9);
  for (auto i = 0; i < F.back().cols(); ++i) {
    Eigen::Vector3i p = F.back().col(i);
    Eigen::VectorXd ints(9);
    ints.setZero();
    // determine all integrals on the element in question
    for (auto j = 0; j < Q._w.size(); ++j) {
      Eigen::Vector2d xi = Q._xi.col(j);
      double w = Q._w(j) * det_dF_T(xi, P.col(p(0)), P.col(p(1)), P.col(p(2)));
      for (auto k = 0; k < 9; ++k)
        ints(k) += w * phi.f(k / 3, xi) * phi.f(k % 3, xi);
    }
    // write integrals to particular entries
    for (auto j = 0; j < 3; ++j)
      for (auto k = 0; k < 3; ++k)
        tripletList.push_back(
            Eigen::Triplet<double>(p(j), p(k), ints(j + 3 * k)));
  }
  M.setFromTriplets(tripletList.begin(), tripletList.end());
  return M;
}

/*
 *  compute all local element matrices on all refinement levels.
 *  finest level ist computed explicitly, coarser levels are obtained
 *  by restriction of the weights on the finer grid.
 */
std::vector<Eigen::MatrixXd> computeMassWeights(const Mesh &msh) {
  auto P = msh.PC();
  auto F = msh.Ffull();
  Phi phi;
  Quadrature Q;
  std::vector<Eigen::MatrixXd> massWeights(F.size());
  Eigen::MatrixXd wght;
  /* restriction weights for coarser levels */
  constexpr double a0[12] = {1., .5, .5, .5, 0., 0., .5, 0., 0., .5, 0., .5};
  constexpr double a1[12] = {0., .5, 0., .5, 1., .5, 0., .5, 0., .5, .5, 0.};
  constexpr double a2[12] = {0., 0., .5, 0., 0., .5, .5, .5, 1., 0., .5, .5};
  auto itsW = massWeights.rbegin();
  for (auto it = F.rbegin(); it != F.rend(); ++it) {
    wght.resize(9, (*it).cols());
    wght.setZero();
    if (it == F.rbegin()) {
      for (auto i = 0; i < it->cols(); ++i)
        wght.col(i) = intMass(P, it->col(i), phi, Q);
    } else {
      for (auto i = 0; i < it->cols(); ++i) {
        auto i4 = 4 * i;
        for (auto j = 0; j < 4; ++j) {
          auto j3 = 3 * j;
          for (auto k = 0; k < 9; ++k) {
            wght(0, i) += a0[j3 + k / 3] * a0[j3 + k % 3] * (*itsW)(k, i4 + j);
            wght(1, i) += a0[j3 + k / 3] * a1[j3 + k % 3] * (*itsW)(k, i4 + j);
            wght(2, i) += a0[j3 + k / 3] * a2[j3 + k % 3] * (*itsW)(k, i4 + j);
            wght(3, i) += a1[j3 + k / 3] * a0[j3 + k % 3] * (*itsW)(k, i4 + j);
            wght(4, i) += a1[j3 + k / 3] * a1[j3 + k % 3] * (*itsW)(k, i4 + j);
            wght(5, i) += a1[j3 + k / 3] * a2[j3 + k % 3] * (*itsW)(k, i4 + j);
            wght(6, i) += a2[j3 + k / 3] * a0[j3 + k % 3] * (*itsW)(k, i4 + j);
            wght(7, i) += a2[j3 + k / 3] * a1[j3 + k % 3] * (*itsW)(k, i4 + j);
            wght(8, i) += a2[j3 + k / 3] * a2[j3 + k % 3] * (*itsW)(k, i4 + j);
          }
        }
      }
      ++itsW;
    }
    *itsW = std::move(wght);
  }
  return massWeights;
}
}  // namespace FEM
}  // namespace SFEMPP
#endif
