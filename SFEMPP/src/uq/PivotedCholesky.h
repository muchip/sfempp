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
#ifndef SFEMPP_UQ_PIVOTEDCHOLESKY_H_
#define SFEMPP_UQ_PIVOTEDCHOLESKY_H_

namespace SFEMPP {
class PivotedCholesky {
 public:
  // void constructor
  PivotedCholesky() {}
  // non-void constructor
  template <class CovKernel>
  PivotedCholesky(const CovKernel &C, double tol) {
    compute(C, tol);
  }
  /*
   *   \brief approximates the approximation error of the pivoted Cholesky
   *   decomposition by Monte Carlo sampling
   */
  template <class CovKernel>
  void sampleError(const CovKernel &C, int samples) const {
    double colOp = 0;
    double colL = 0;
    double maxError = 0;
    double aveError = 0;
    double error = 0;
    int col = 0;
    int dim = C.dim();
    std::srand(std::time(NULL));
    // sample random columns from C and compare them to the respective column
    // from L * L'
    for (int i = 0; i < samples; ++i) {
      error = 0;
      col = std::rand() % dim;
      for (int j = 0; j < dim; ++j) {
        colOp = C(j, col);
        colL = L_.row(j).dot(L_.row(col));
        error += (colOp - colL) * (colOp - colL);
      }
      error = sqrt(error);
      if (maxError < error) maxError = error;
      aveError += error / samples;
    }

    std::cout << "maximum error: " << maxError << " average error: " << aveError
              << " samples: " << samples << std::endl;
  }
  const Eigen::MatrixXd &matrixL(void) const { return L_; }
  double tolerance(void) const { return tol_; }
  /*
   *    \brief computes the pivoted Cholesky decomposition of a positive
   *           semi definite function provided by the class covKernel
   */
  template <class CovKernel>
  void compute(const CovKernel &C, double tol) {
    Eigen::VectorXd D;
    double tr = 0;
    constexpr int allocBSize = 100;
    int actBSize = 0;
    int step = 0;
    int pivot = 0;
    int dim = C.dim();
    tol_ = tol;
    D.resize(dim);
    // compute diagonal
    for (int i = 0; i < dim; ++i) D(i) = C(i, i);
    // scaling prec by the trace for relative error
    tr = D.sum();
    tol *= tr;
    // perform pivoted Cholesky decomposition
    for (step = 0; (step < dim) && (tol < tr); ++step) {
      D.maxCoeff(&pivot);
      if (actBSize - 1 <= step) {
        actBSize += allocBSize;
        L_.conservativeResize(dim, actBSize);
      }
      L_.col(step).setZero();
      for (int i = 0; i < dim; ++i) L_(i, step) = C(i, pivot);
      L_.col(step) -= L_.leftCols(step) * L_.row(pivot).head(step).transpose();
      if (L_(pivot, step) <= 0) {
        std::cout << "breaking with non positive pivot\n";
        break;
      }
      L_.col(step) /= sqrt(L_(pivot, step));
      // update the diagonal and the trace
      D.array() -= L_.col(step).array().square();
      // update the trace
      tr = D.sum();
      ++step;
    }
    // crop L to its actual size
    L_.conservativeResize(dim, step);
  }

 private:
  double tol_;
  Eigen::MatrixXd L_;
};
}  // namespace SFEMPP

#endif
