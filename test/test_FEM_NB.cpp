#include "../SFEMPP/FEM"
#include "../SFEMPP/MMG"
#include "../SFEMPP/UQ"
#include "../SFEMPP/src/util/TicToc.h"

struct func {
  double operator()(const Eigen::Vector2d &x) const {
    return 100 * x.norm() * exp(-10 * x.norm());
  }
};

int main(int argc, char *argv[]) {
  SFEMPP::TicToc T;
  int maxlvl = atoi(argv[1]);
  // geometry and mesh
  T.tic();
  SFEMPP::Disc domain;
  SFEMPP::Mesh myMesh(domain);
  myMesh.refine(maxlvl).setCurved(domain);
  myMesh.setNeumann();
  T.toc("mesh generation:                   ");
  std::cout << "number of points:                   " << myMesh.P().cols()
            << std::endl;
  std::cout << "number of elements:                 " << myMesh.F(maxlvl).cols()
            << std::endl;
  T.tic();
  SFEMPP::CovarianceKernelMF cov_kernel(myMesh.PC());
  SFEMPP::KarhunenLoeveMF kl;
  kl.compute(cov_kernel, myMesh, 1e-5);
  T.toc("computation of the KL expansion: ");
  kl.prune(1e-5);
  kl.sampleError(cov_kernel, 10);
  std::cout << "size kl: " << kl.KL().cols() << std::endl;
  Eigen::MatrixXd PCref = myMesh.PC();
  Eigen::MatrixXd V = 0.3 * kl.KL() * Eigen::VectorXd::Random(kl.KL().cols());
  myMesh.set_PC(PCref + V.reshaped(PCref.cols(), 2).transpose());
  std::vector<Eigen::MatrixXd> stiffness_weights =
      SFEMPP::FEM::computeStiffnessWeights(myMesh);
  std::vector<Eigen::MatrixXd> mass_weights =
      SFEMPP::FEM::computeMassWeights(myMesh);
  T.toc("computation of stiffness weights:  ");
  Eigen::VectorXd d = Eigen::VectorXd::Ones(myMesh.F(maxlvl).cols());
  T.tic();
  auto S = SFEMPP::FEM::computeStiffnessMatrixMG(myMesh, stiffness_weights, d);
  auto M = SFEMPP::FEM::computeMassMatrixMG(myMesh, mass_weights);
  auto SpM = S;
  for (int i = 0; i < M.size(); ++i) SpM[i] += M[i];
  T.toc("computation of a stiffness matrix: ");
  auto f = SFEMPP::FEM::computeLoadVector(myMesh, func());
  Eigen::VectorXd uref(myMesh.nP()[maxlvl]);
  uref.setZero();
  T.tic();
  SFEMPP::MMG::mmg(SpM, f, uref, myMesh, maxlvl);
  myMesh.printSolution(uref, "sol.vtk");
  //////////////////////////////////////////////////////////////////////////////
  std::vector<Eigen::VectorXd> rhs(maxlvl + 1);
  rhs.back() = f;
  for (int l = maxlvl - 1; l >= 0; --l) {
    rhs[l] = rhs[l + 1];
    SFEMPP::MMG::restriction(rhs[l], myMesh, l + 1);
  }
  for (auto l = 0; l < maxlvl; ++l) {
    T.tic();
    Eigen::VectorXd ul(myMesh.nP()[l]);
    Eigen::VectorXd dl = Eigen::VectorXd::Ones(myMesh.F(l).cols());
    ul.setZero();
    auto Sl =
        SFEMPP::FEM::computeStiffnessMatrixMG(myMesh, stiffness_weights, d, l);
    for (int i = 0; i < Sl.size(); ++i) Sl[i] += M[i];

    SFEMPP::MMG::mmg(Sl, rhs[l], ul, myMesh, l);
    ul = SFEMPP::MMG::projection(myMesh, ul, l, maxlvl);
    T.toc("solution on level " + std::to_string(l) + ":         ");
    std::cout << "L2 error: "
              << SFEMPP::FEM::computeL2norm(myMesh.PC(), myMesh.F(maxlvl),
                                            ul - uref)
              << std::endl;
  }

  return 0;
}
