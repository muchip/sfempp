#include "../SFEMPP/FEM"
#include "../SFEMPP/MMG"
#include "../SFEMPP/src/util/TicToc.h"

struct func {
  double operator()(const Eigen::Vector2d &x) const { return 4.; }
};

struct uref {
  double operator()(const Eigen::Vector2d &x) const {
    return 1 - x(0) * x(0) - x(1) * x(1);
  }
};

struct grad_uref {
  Eigen::Vector2d operator()(const Eigen::Vector2d &x) const {
    Eigen::Vector2d retval;
    retval << -2 * x(0), -2 * x(1);
    return retval;
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
  T.toc("mesh generation:                   ");
  std::cout << "number of points:                   " << myMesh.P().cols()
            << std::endl;
  std::cout << "number of elements:                 " << myMesh.F(maxlvl).cols()
            << std::endl;
  T.tic();
  std::vector<Eigen::MatrixXd> stiffness_weights =
      SFEMPP::FEM::computeStiffnessWeights(myMesh);
  T.toc("computation of stiffness weights:  ");
  Eigen::VectorXd d = Eigen::VectorXd::Ones(myMesh.F(maxlvl).cols());
  T.tic();
  auto S = SFEMPP::FEM::computeStiffnessMatrixMG(myMesh, stiffness_weights, d);
  T.toc("computation of a stiffness matrix: ");
  auto f = SFEMPP::FEM::computeLoadVector(myMesh, func());
  Eigen::VectorXd u(myMesh.nP()[maxlvl]);
  u.setZero();
  T.tic();
  SFEMPP::MMG::mmg(S, f, u, myMesh, maxlvl);
  T.toc("computation of the solution:       ");
  double error = SFEMPP::FEM::computeL2error(myMesh, u, uref());
  std::cout << "L2-error:                           " << error << std::endl;
  error = SFEMPP::FEM::computeH1error(myMesh, u, grad_uref());
  std::cout << "H1-error:                           " << error << std::endl;
  myMesh.printSolution(u, "sol.vtk");

  return 0;
}
