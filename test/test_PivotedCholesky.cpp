
#include "../SFEMPP/Mesh"
#include "../SFEMPP/UQ"
#include "../SFEMPP/src/util/TicToc.h"

int main(int argc, char *argv[]) {
  SFEMPP::TicToc T;
  int maxlvl = 9;
  // geometry and mesh
  T.tic();
  SFEMPP::Kite domain;
  SFEMPP::Mesh myMesh(domain);
  myMesh.refine(maxlvl).setCurved(domain);
  T.toc("computation of the geometry:   ");
  std::cout << "number of points:               " << myMesh.P().cols()
            << std::endl;
  std::cout << "number of elements:             " << myMesh.F(maxlvl).cols()
            << std::endl;
  Eigen::Vector2d minP = myMesh.P().rowwise().minCoeff();
  Eigen::Vector2d maxP = myMesh.P().rowwise().maxCoeff();
  Eigen::Vector2d diamP = maxP - minP;
  std::cout << "bounding box size:              " << diamP(0) << " x "
            << diamP(1) << std::endl;
  SFEMPP::CovarianceKernelMF cov_kernel(myMesh.P());
  std::cout << "kernel matrix size:             " << cov_kernel.dim() << " x "
            << cov_kernel.dim() << std::endl;
  T.tic();
  SFEMPP::PivotedCholesky pc;
  pc.compute(cov_kernel, 1e-10);
  T.toc("computation of the Pchol:      ");
  std::cout << "rank of Pchol:                  " << pc.matrixL().cols()
            << std::endl;
  pc.sampleError(cov_kernel, 10);
  return 0;
}
