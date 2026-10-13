#include "../SFEMPP/FEM"
#include "../SFEMPP/MMG"
#include "../SFEMPP/UQ"
#include "../SFEMPP/src/util/TicToc.h"

struct func {
  double operator()(const Eigen::Vector2d &x) const {
    return 1. / (std::sqrt((x(0) - 1) * (x(0) - 1) + x(1) * x(1)) + 1e-16);
  }
};

int main(int argc, char *argv[]) {
  SFEMPP::TicToc T;
  int maxlvl = atoi(argv[1]);
  // geometry and mesh
  T.tic();
  SFEMPP::Kite domain;
  SFEMPP::Mesh myMesh(domain);
  myMesh.refine(maxlvl).setCurved(domain);
  myMesh.setNeumann();
  T.toc("mesh generation:                   ");
  std::cout << "number of points:                   " << myMesh.P().cols()
            << std::endl;
  std::cout << "number of elements:                 " << myMesh.F(maxlvl).cols()
            << std::endl;
  Eigen::MatrixXd PCref = myMesh.PC();
  PCref *= 0.2;
  myMesh.set_PC(PCref);

  T.tic();
  SFEMPP::CovarianceKernelMF cov_kernel(myMesh.PC());
  SFEMPP::KarhunenLoeveMF kl;
  kl.compute(cov_kernel, myMesh, 1e-5);
  T.toc("computation of the KL expansion: ");
  kl.prune(1e-5);
  kl.sampleError(cov_kernel, 10);
  std::cout << "size kl: " << kl.KL().cols() << std::endl;
  Eigen::VectorXd uref(myMesh.nP()[maxlvl]);
  T.tic();
  Eigen::VectorXd mean = uref;
  Eigen::VectorXd var = uref;
  mean.setZero();
  var.setZero();
  for (int sample = 0; sample < 1e4; ++sample) {
    Eigen::MatrixXd V = 0.2 * kl.KL() * Eigen::VectorXd::Random(kl.KL().cols());
    myMesh.set_PC(PCref + V.reshaped(PCref.cols(), 2).transpose());
    std::vector<Eigen::MatrixXd> stiffness_weights =
        SFEMPP::FEM::computeStiffnessWeights(myMesh);
    std::vector<Eigen::MatrixXd> mass_weights =
        SFEMPP::FEM::computeMassWeights(myMesh);
    Eigen::VectorXd d = Eigen::VectorXd::Ones(myMesh.F(maxlvl).cols());
    auto S =
        SFEMPP::FEM::computeStiffnessMatrixMG(myMesh, stiffness_weights, d);
    auto M = SFEMPP::FEM::computeMassMatrixMG(myMesh, mass_weights);
    auto SpM = S;
    for (int i = 0; i < M.size(); ++i) SpM[i] += M[i];
    auto f = SFEMPP::FEM::computeLoadVector(myMesh, func());
    uref.setZero();
    SFEMPP::MMG::mmg(SpM, f, uref, myMesh, maxlvl);
    mean += 1e-4 * uref;
    var += (1e-4 * uref.array().square()).matrix();
    if (sample < 50)
      myMesh.printSolution(uref, "sol" + std::to_string(sample) + ".vtk");
  }
  T.toc("sampling: ");
  myMesh.set_PC(PCref);
  myMesh.printSolution(mean, "mean.vtk");
  myMesh.printSolution(var.array() - mean.array().square(), "var.vtk");

  return 0;
}
