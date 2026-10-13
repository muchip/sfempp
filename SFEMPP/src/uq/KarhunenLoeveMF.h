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
#ifndef SFEMPP_UQ_KARHUNENLOEVEMF_H_
#define SFEMPP_UQ_KARHUNENLOEVEMF_H_

namespace SFEMPP {
class KarhunenLoeveMF {
 public:
  KarhunenLoeveMF() {}
  template <class covKernel>
  KarhunenLoeveMF(const covKernel &C, const Mesh &myMesh, double tol) {
    compute(C, myMesh, tol);
  }

  const Eigen::MatrixXd &KL() const { return KL_; }

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
        colL = KL_.row(j).dot(KL_.row(col));
        error += (colOp - colL) * (colOp - colL);
      }
      error = sqrt(error);
      if (maxError < error) maxError = error;
      aveError += error / samples;
    }

    std::cout << "maximum error: " << maxError << " average error: " << aveError
              << " samples: " << samples << std::endl;
  }

  template <class covKernel>
  void compute(const covKernel &C, const Mesh &myMesh, double tol) {
    PivotedCholesky pC(C, tol * tol);
    unsigned int n = myMesh.P().cols();
    auto M = FEM::computeMassMatrix(myMesh);
    Eigen::MatrixXd MML(pC.matrixL().rows(), pC.matrixL().cols());
    MML.topRows(n) = M * pC.matrixL().topRows(n);
    MML.bottomRows(n) = M * pC.matrixL().bottomRows(n);
    Eigen::MatrixXd Cmat = pC.matrixL().transpose() * MML;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Cmat);
    KL_ = pC.matrixL() * es.eigenvectors();
    for (auto i = 0; i < KL_.cols() / 2; ++i)
      KL_.col(i).swap(KL_.col(KL_.cols() - 1 - i));

    gamma_.resize(KL_.cols());
    for (auto i = 0; i < gamma_.size(); ++i)
      gamma_(i) = KL_.col(i).cwiseAbs().maxCoeff();
  }
  const Eigen::VectorXd &gamma() const { return gamma_; }

  void prune(double tol) {
    Eigen::MatrixXd KL;
    Eigen::VectorXi sortVec =
        Eigen::ArrayXi::LinSpaced(gamma_.size(), 0, gamma_.size() - 1);
    std::sort(sortVec.data(), sortVec.data() + sortVec.size(),
              myComp<Eigen::VectorXd>(gamma_));
    double sumGamma = gamma_.sum();
    double curSum = 0;
    int numVec = 0;
    while (sumGamma - curSum > tol * sumGamma) {
      curSum += gamma_(sortVec(numVec));
      ++numVec;
    }
    KL.resize(KL_.rows(), numVec);
    for (auto i = 0; i < numVec; ++i) KL.col(i) = KL_.col(sortVec(i));
    KL_ = KL;
    gamma_.resize(KL_.cols());
    for (auto i = 0; i < gamma_.size(); ++i)
      gamma_(i) = KL_.col(i).cwiseAbs().maxCoeff();
  }

 protected:
  // comparison functors for C++ sort routine
  template <typename T>
  struct myComp {
    const T &m;
    myComp(T &p) : m(p){};

    bool operator()(const int &i, const int &j) { return (m(i) > m(j)); }
  };
  Eigen::MatrixXd KL_;
  Eigen::VectorXd gamma_;
};
}  // namespace SFEMPP
#endif
