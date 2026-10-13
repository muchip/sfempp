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

#ifndef SFEMPP_FEM_STIFFNESSMATRIX_H_
#define SFEMPP_FEM_STIFFNESSMATRIX_H_

namespace SFEMPP {
namespace FEM {
/*
 *  computes all integrals a(phi_i,phi_j) on the patch
 * (P[F[0]],P[F[1]],P[F[2]])
 *  we use only the polygonal approxmation here and fix the quadrature to the
 *  radon rule
 */
// constant coefficient version
Eigen::VectorXd intStiff(const Eigen::MatrixXd &P, const Eigen::VectorXi &F,
                         const Phi &phi, const Quadrature &Q) {
  using namespace Transformations;
  Eigen::VectorXd stiff(9);
  stiff.setZero();
  for (auto i = 0; i < Q._w.size(); ++i) {
    Eigen::Vector2d xi = Q._xi.col(i);
    Eigen::Vector2d a = dF_Tdx(xi, P.col(F(0)), P.col(F(1)), P.col(F(2)));
    Eigen::Vector2d b = dF_Tdy(xi, P.col(F(0)), P.col(F(1)), P.col(F(2)));
    double c = det_dF_T(xi, P.col(F(0)), P.col(F(1)), P.col(F(2)));
    Eigen::Vector2d dPhi1;
    Eigen::Vector2d dPhi2;
    for (auto k = 0; k < 9; ++k) {
      dPhi1(0) = phi.dx(k / 3, xi) * b(1) - phi.dy(k / 3, xi) * a(1);
      dPhi1(1) = phi.dy(k / 3, xi) * a(0) - phi.dx(k / 3, xi) * b(0);
      dPhi2(0) = phi.dx(k % 3, xi) * b(1) - phi.dy(k % 3, xi) * a(1);
      dPhi2(1) = phi.dy(k % 3, xi) * a(0) - phi.dx(k % 3, xi) * b(0);
      // note that the term /c comes from Cramer's rule for the transported
      // gradient onto the reference element.
      stiff(k) += Q._w(i) * (dPhi1(0) * dPhi2(0) + dPhi1(1) * dPhi2(1)) / c;
    }
  }
  return stiff;
}
// variable coefficient version
Eigen::VectorXd intStiff(const Eigen::Matrix2d &A, const Eigen::MatrixXd &P,
                         const Eigen::VectorXi &F, const Phi &phi,
                         const Quadrature &Q) {
  using namespace Transformations;
  Eigen::VectorXd stiff(9);
  stiff.setZero();
  for (auto i = 0; i < Q._w.size(); ++i) {
    Eigen::Vector2d xi = Q._xi.col(i);
    Eigen::Vector2d a = dF_Tdx(xi, P.col(F(0)), P.col(F(1)), P.col(F(2)));
    Eigen::Vector2d b = dF_Tdy(xi, P.col(F(0)), P.col(F(1)), P.col(F(2)));
    double c = det_dF_T(xi, P.col(F(0)), P.col(F(1)), P.col(F(2)));
    Eigen::Vector2d dPhi1;
    Eigen::Vector2d dPhi2;
    for (auto k = 0; k < 9; ++k) {
      dPhi1(0) = phi.dx(k / 3, xi) * b(1) - phi.dy(k / 3, xi) * a(1);
      dPhi1(1) = phi.dy(k / 3, xi) * a(0) - phi.dx(k / 3, xi) * b(0);
      dPhi2(0) = phi.dx(k % 3, xi) * b(1) - phi.dy(k % 3, xi) * a(1);
      dPhi2(1) = phi.dy(k % 3, xi) * a(0) - phi.dx(k % 3, xi) * b(0);
      // note that the term /c comes from Cramer's rule for the transported
      // gradient onto the reference element.
      stiff(k) += Q._w(i) * (A * dPhi1).dot(dPhi2) / c;
    }
  }
  return stiff;
}
/*
 *  compute all local element matrices on all refinement levels.
 *  finest level ist computed explicitly, coarser levels are obtained
 *  by restriction of the weights on the finer grid.
 */
std::vector<Eigen::MatrixXd> computeStiffnessWeights(const Mesh &msh) {
  auto P = msh.PC();
  auto F = msh.Ffull();
  Phi phi;
  Quadrature Q;
  std::vector<Eigen::MatrixXd> stiffnessWeights(F.size());
  Eigen::MatrixXd wght;
  /* restriction weights for coarser levels */
  constexpr double a0[12] = {1., .5, .5, .5, 0., 0., .5, 0., 0., .5, 0., .5};
  constexpr double a1[12] = {0., .5, 0., .5, 1., .5, 0., .5, 0., .5, .5, 0.};
  constexpr double a2[12] = {0., 0., .5, 0., 0., .5, .5, .5, 1., 0., .5, .5};
  auto itsW = stiffnessWeights.rbegin();
  for (auto it = F.rbegin(); it != F.rend(); ++it) {
    wght.resize(9, (*it).cols());
    wght.setZero();
    if (it == F.rbegin()) {
      for (auto i = 0; i < it->cols(); ++i)
        wght.col(i) = intStiff(P, it->col(i), phi, Q);
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
  return stiffnessWeights;
}

std::vector<Eigen::SparseMatrix<double>> computeStiffnessMatrixMG(
    const Mesh &msh, const std::vector<Eigen::MatrixXd> &stiffnessWeights,
    const Eigen::VectorXd &d, int l = -1) {
  auto maxLvl = msh.maxlvl();
  auto P = msh.PC();
  auto nP = msh.nP();
  auto F = msh.Ffull();
  auto B = msh.B();
  Eigen::VectorXd myD = d;
  Eigen::Vector3i p;
  if (l == -1 || l > maxLvl) l = maxLvl;
  std::vector<Eigen::SparseMatrix<double>> stiffnessMatrix(l + 1);
  auto itA = stiffnessMatrix.rbegin();
  auto itF = F.rbegin();
  auto itn = nP.rbegin();
  auto itW = stiffnessWeights.rbegin();
  std::advance(itF, maxLvl - l);
  std::advance(itn, maxLvl - l);
  std::advance(itW, maxLvl - l);
  for (; itF != F.rend(); ++itA, ++itF, ++itn, ++itW) {
    Eigen::SparseMatrix<double> S(*itn, *itn);
    S.reserve((*itn) * (msh.computeMaxDegree() + 1));
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve((*itn) * (msh.computeMaxDegree() + 1) * 9);
    // add entries into sparse matrix
    for (auto i = 0; i < itF->cols(); ++i) {
      p = itF->col(i);
      for (auto j = 0; j < 3; ++j)
        for (auto k = 0; k < 3; ++k)
          if (B(p(j)) && B(p(k)))
            tripletList.push_back(Eigen::Triplet<double>(
                p(j), p(k), myD(i) * (*itW)(j + 3 * k, i)));
    }
    // restrict diffusion coefficient to next level
    for (auto i = 0; i < itF->cols() / 4; ++i)
      myD(i) = 0.25 * myD.segment(4 * i, 4).sum();
    S.setFromTriplets(tripletList.begin(), tripletList.end());
    *itA = std::move(S);
  }
  return stiffnessMatrix;
}

/**
 *   \brief computes the stiffness matrix with an anisotropic coefficient
 **/
Eigen::SparseMatrix<double> computeStiffnessMatrix(const Mesh &msh,
                                                   const Eigen::MatrixXd &D) {
  auto P = msh.PC();
  auto E = msh.Ffull().back();
  auto B = msh.B();
  Eigen::Vector3i p;
  Eigen::Matrix2d A;
  Phi phi;
  Quadrature Q;
  Eigen::SparseMatrix<double> S(P.cols(), P.cols());
  S.reserve(P.cols() * (msh.computeMaxDegree() + 1));
  std::vector<Eigen::Triplet<double>> tripletList;
  tripletList.reserve(P.cols() * (msh.computeMaxDegree() + 1) * 9);
  // add entries into sparse matrix
  for (auto i = 0; i < E.cols(); ++i) {
    p = E.col(i);
    A << D(0, i), D(1, i), D(2, i), D(3, i);
    auto ints = intStiff(A, P, p, phi, Q);
    for (auto j = 0; j < 3; ++j)
      for (auto k = 0; k < 3; ++k)
        if (B(p(j)) && B(p(k)))
          tripletList.push_back(
              Eigen::Triplet<double>(p(j), p(k), ints(j + 3 * k)));
  }
  S.setFromTriplets(tripletList.begin(), tripletList.end());
  return S;
}
}  // namespace FEM
}  // namespace SFEMPP
#endif
