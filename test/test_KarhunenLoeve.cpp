
#include "../SFEMPP/Mesh"
#include "../SFEMPP/UQ"
#include "../SFEMPP/src/util/TicToc.h"

int main(int argc, char *argv[]) {
  SFEMPP::TicToc T;
  int maxlvl = 7;
  // geometry and mesh
  T.tic();
  SFEMPP::Star domain;
  SFEMPP::Mesh myMesh(domain);
  myMesh.refine(maxlvl).setCurved(domain);
  T.toc("computation of the geometry:   ");
  std::cout << "number of points:               " << myMesh.P().cols()
            << std::endl;
  std::cout << "number of elements:             " << myMesh.F(maxlvl).cols()
            << std::endl;
  Eigen::Vector2d minP = myMesh.PC().rowwise().minCoeff();
  Eigen::Vector2d maxP = myMesh.PC().rowwise().maxCoeff();
  Eigen::Vector2d diamP = maxP - minP;
  std::cout << "bounding box size:              " << diamP(0) << " x "
            << diamP(1) << std::endl;
  SFEMPP::CovarianceKernelMF cov_kernel(myMesh.PC());
  std::cout << "kernel matrix size:             " << cov_kernel.dim() << " x "
            << cov_kernel.dim() << std::endl;
  T.tic();
  SFEMPP::KarhunenLoeveMF kl;
  kl.compute(cov_kernel, myMesh, 1e-5);
  T.toc("computation of the KL expansion: ");
  kl.sampleError(cov_kernel, 10);
  std::cout << "size kl: " << kl.KL().cols() << std::endl;
  kl.prune(1e-5);
  std::cout << "size kl after pruning: " << kl.KL().cols() << std::endl;
  kl.sampleError(cov_kernel, 10);
  Eigen::MatrixXd PCref = myMesh.PC();
  for (int i = 0; i < 10; ++i) {
    Eigen::MatrixXd V = 0.3 * kl.KL() * Eigen::VectorXd::Random(kl.KL().cols());
    myMesh.set_PC(PCref + V.reshaped(PCref.cols(), 2).transpose());
    myMesh.printGeometry(Eigen::VectorXd::Zero(myMesh.nP().back()), maxlvl,
                         "Kite" + std::to_string(i) + ".vtk");
  }
  return 0;
}
