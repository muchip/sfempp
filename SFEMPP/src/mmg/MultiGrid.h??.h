#ifndef __MULTIGRID__
#define __MULTIGRID__

#include <iomanip>
#include <iostream>
#include "../../Mesh.h"
#include "../Sparse.h"

namespace MG {

const int __minLvl = 1;
const int __cycleType = 2;
const int __preSmooth = 3;
const int __postSmooth = 3;
/*
 *  computs the restriction of the vector f from level m to m-1
 */
void restriction(Eigen::VectorXd &f, const Mesh &msh, int m) {
  const Eigen::MatrixXi &F = msh.get_F(m);
  const Eigen::MatrixXi &resF = msh.get_F(m - 1);
  Eigen::VectorXi p;
  Eigen::VectorXi mp;
  Eigen::Vector3d a;
  for (auto i = 0; i < resF.cols(); ++i) {
    // determine vertices of the element on level m-1
    p = resF.col(i);
    mp = F.col(4 * i + 3);
    // determine values at the midpoints
    for (auto j = 0; j < 3; ++j) a(j) = msh.get_B()(mp(j)) * f(mp(j));
    // restriction
    for (auto j = 0; j < 3; ++j)
      if (msh.get_B()(p(j))) f(p(j)) += 0.25 * (a(j) + a((j + 2) % 3));
  }
  f.conservativeResize(msh.get_nP()[m - 1]);
}

/*
 *  computs the prolongation of the vector f from level m-1 to m
 */
void prolongation(Eigen::VectorXd &f, const Mesh &msh, int m) {
  const Eigen::MatrixXi &F = msh.get_F(m - 1);
  const Eigen::MatrixXi &proF = msh.get_F(m);
  Eigen::VectorXi mp;
  Eigen::Vector3d a;
  f.conservativeResize(msh.get_nP()[m]);
  f.tail(msh.get_nP()[m] - msh.get_nP()[m - 1]).setZero();
  for (auto i = 0; i < F.cols(); ++i) {
    // get values at the vertices of the element
    for (auto j = 0; j < 3; ++j) a(j) = f(F(j, i));
    mp = proF.col(4 * i + 3);
    // distribute values to the midpoints
    for (auto j = 0; j < 3; ++j)
      f(mp(j)) += 0.25 * msh.get_B()(mp(j)) * (a(j) + a((j + 1) % 3));
  }
}
/*
 *  Add Jacobi smoother for use with multigrid. Gauss Seidel will probably
 *  not work efficiently...
 */
void JacobiSmoother(const Sparse<double> &M,
                    const Eigen::VectorXd &b, Eigen::VectorXd &x, int steps) {
  Eigen::VectorXd D = M.diagonal();
  for (auto i = 0; i < steps; ++i) {
    // multiply with full sparse matrix
    x = b - (M * x - x.cwiseProduct(D));
    for (auto j = 0; j < x.size(); ++j)
      if (D(j))
        x(j) /= D(j);
      else
        x(j) = 0;
    // compute difference, repair diagonal and add scaling
    // x = (b.array() - x.array() + D.array() * x.array()) / D.array();
    for (auto j = 0; j < x.size(); ++j)
      if (x(j) != x(j)) x(j) = 0;
  }
}

/*
 *  Add Jacobi smoother for use with multigrid. Gauss Seidel will probably
 *  not work efficiently...
 */
void GaussSeidelSmoother(const Sparse<double> &M,
                         const Eigen::VectorXd &b, Eigen::VectorXd &x,
                         int steps) {
  double c = 0;
  for (auto i = 0; i < steps; ++i)
    for (auto j = 0; j < x.size(); ++j) {
      x(j) = 0;
      for (auto k = 0; k < M.row(j).size(); ++k)
        if (M.row(j)[k].index != j)
          x(j) += M.row(j)[k].value * x(M.row(j)[k].index);
        else
          c = M.row(j)[k].value;
      x(j) = (b(j) - x(j)) / c;
    }
}

void mmgStep(const std::vector<Sparse<double> > &M,
             const Eigen::VectorXd &f, Eigen::VectorXd &x, const Mesh &msh,
             int lvl) {
  if (lvl <= __minLvl)
    x = M[lvl].full().fullPivLu().solve(f);
  else {
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

void mmg(const std::vector<Sparse<double> > &M,
         const Eigen::VectorXd &f, Eigen::VectorXd &x, const Mesh &msh,
         int lvl) {
  Eigen::VectorXd res = f - M[lvl] * x;
  int i = 0;
  while (res.norm() > 1e-6) {
    mmgStep(M, f, x, msh, lvl);
    res = f - M[lvl] * x;
    ++i;
  }
  std::cout << "MG steps: " << i << std::endl;
}
}  // namespace MG
#endif
