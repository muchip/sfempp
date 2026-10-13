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

#ifndef SFEMPP_MMG_MULTIGRID_H_
#define SFEMPP_MMG_MULTIGRID_H_

namespace SFEMPP {
namespace MMG {

const int __minLvl = 1;
const int __cycleType = 1;
const int __preSmooth = 3;
const int __postSmooth = 3;

/*
 *  computes the restriction of the vector f from level m to m-1
 */
void restriction(Eigen::VectorXd &f, const Mesh &msh, int m) {
  auto B = msh.B();
  auto F = msh.F(m);
  auto resF = msh.F(m - 1);
  Eigen::VectorXi p;
  Eigen::VectorXi mp;
  Eigen::Vector3d a;
  for (auto i = 0; i < resF.cols(); ++i) {
    // determine vertices of the element on level m-1
    p = resF.col(i);
    mp = F.col(4 * i + 3);
    // determine values at the midpoints
    for (auto j = 0; j < 3; ++j) a(j) = B(mp(j)) * f(mp(j));
    // restriction
    for (auto j = 0; j < 3; ++j)
      if (B(p(j))) f(p(j)) += 0.25 * (a(j) + a((j + 2) % 3));
  }
  f.conservativeResize(msh.nP()[m - 1]);
}

/*
 *  computes the prolongation of the vector f from level m-1 to m
 */
void prolongation(Eigen::VectorXd &f, const Mesh &msh, int m) {
  auto B = msh.B();
  auto F = msh.F(m - 1);
  auto proF = msh.F(m);
  Eigen::VectorXi mp;
  Eigen::Vector3d a;
  f.conservativeResize(msh.nP()[m]);
  f.tail(msh.nP()[m] - msh.nP()[m - 1]).setZero();
  for (auto i = 0; i < F.cols(); ++i) {
    // get values at the vertices of the element
    for (auto j = 0; j < 3; ++j) a(j) = f(F(j, i));
    mp = proF.col(4 * i + 3);
    // distribute values to the midpoints
    for (auto j = 0; j < 3; ++j)
      f(mp(j)) += 0.25 * B(mp(j)) * (a(j) + a((j + 1) % 3));
  }
}

// buggy since there are zeros on the diagonal of M!!!
// however, it works as long as the right hand side is exactly zero in
// those rows due to the sparse optimization in the eigen solve method
void GaussSeidelSmoother(const Eigen::SparseMatrix<double> &M,
                         const Eigen::VectorXd &b, Eigen::VectorXd &x,
                         int steps) {
  for (auto i = 0; i < steps; ++i)
    x = M.triangularView<Eigen::Lower>().solve(
        b - M.triangularView<Eigen::StrictlyUpper>() * x);
}

void mmgStep(const std::vector<Eigen::SparseMatrix<double> > &M,
             const Eigen::VectorXd &f, Eigen::VectorXd &x, const Mesh &msh,
             int lvl) {
  if (lvl <= __minLvl) {
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >
        solver;
    solver.compute(M[lvl]);
    x = solver.solve(f);
  } else {
    GaussSeidelSmoother(M[lvl], f, x, __preSmooth);
    Eigen::VectorXd r = f - M[lvl] * x;
    restriction(r, msh, lvl);
    Eigen::VectorXd y = Eigen::VectorXd::Zero(r.size());
    for (auto i = 0; i < __cycleType; ++i) mmgStep(M, r, y, msh, lvl - 1);
    prolongation(y, msh, lvl);
    x += y;
    GaussSeidelSmoother(M[lvl], f, x, __postSmooth);
  }
}

void mmg(const std::vector<Eigen::SparseMatrix<double> > &M,
         const Eigen::VectorXd &f, Eigen::VectorXd &x, const Mesh &msh, int lvl,
         double tol = 1e-6) {
  int iter = 0;
  Eigen::VectorXd res = f - M[lvl] * x;
  while (res.norm() > tol) {
    mmgStep(M, f, x, msh, lvl);
    res = f - M[lvl] * x;
    ++iter;
  }
  // std::cout << "MG iter: " << iter << std::endl;
}

}  // namespace MMG
}  // namespace SFEMPP
#endif
