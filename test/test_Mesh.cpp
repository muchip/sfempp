#include "../SFEMPP/Mesh"
#include "../SFEMPP/src/util/TicToc.h"

int main(int argc, char *argv[]) {
  SFEMPP::TicToc T;
  int maxlvl = 6;
  // geometry and mesh
  T.tic();
  SFEMPP::UQshape domain;
  SFEMPP::Mesh myMesh(domain);
  myMesh.refine(maxlvl).setCurved(domain);
  myMesh.printGeometry(Eigen::VectorXd::Zero(myMesh.nP().back()), maxlvl,
                       "Kite.vtk");
  // myMesh.setNeumann();
  return 0;
}
