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
#ifndef SFEMPP_MESH_MESH_H_
#define SFEMPP_MESH_MESH_H_

namespace SFEMPP {
class Mesh {
 public:
  // constructors
  Mesh(const ParametrisationBase &Chi) { genTriang(Chi); }
  // getters
  const Eigen::MatrixXi &F(int m) const { return _F[m]; }
  const std::vector<Eigen::MatrixXi> &Ffull(void) const { return _F; }
  const std::vector<int> &nP(void) const { return _nP; }
  const Eigen::MatrixXd &P(void) const { return _P; }
  const Eigen::MatrixXd &PC(void) const { return _PC; }
  const Eigen::VectorXi &B(void) const { return _B; }
  const unsigned int maxlvl(void) const { return _maxLvl; }
  // setter
  void set_PC(const Eigen::MatrixXd &PC) { _PC = PC; }
  void swap_PC(Eigen::MatrixXd &PC) { _PC.swap(PC); }
  // transform polygonal mesh wrt reference patches to curved
  void setCurved(const ParametrisationBase &Chi) {
    _PC = polygonal2curved(_P, _F, Chi);
  }
  // refine current mesh
  Mesh &refine(void) {
    Eigen::MatrixXi E = getEdgeList(_F.back(), _P);
    // count number of edges to determine memory requirements for new points
    int nEdges = 0;
    for (auto i = 0; i < E.rows(); ++i)
      for (auto j = 0; E(i, j) != -1; ++j) ++nEdges;
    // get memory for new mesh
    Eigen::MatrixXd newP(_P.rows(), nEdges);
    Eigen::VectorXi newB(nEdges);
    Eigen::MatrixXi newF(_F.back().rows(), 4 * _F.back().cols());
    Eigen::MatrixXi newPind = E;
    int currPoint = 0;
    // determine new points
    for (auto i = 0; i < E.rows(); ++i)
      for (auto j = 0; E(i, j) != -1; ++j) {
        // generate new point in as the middle of the current edge
        newP.col(currPoint) = 0.5 * (_P.col(i) + _P.col(E(i, j)));
        // set boundary condition (if both points are at the boundary set
        // the new point also as boundary point)
        if (_B(i) != 1 && _B(E(i, j)) != 1)
          newB(currPoint) = _B(i);
        else
          newB(currPoint) = 1;
        // store location of the new point in a similar pattern as E
        newPind(i, j) = currPoint;
        ++currPoint;
      }
    // update element list
    for (auto i = 0; i < _F.back().cols(); ++i) {
      // determine indices of the midpoints of the current element
      int m[3] = {0, 0, 0};
      for (auto j = 0; j < 3; ++j) {
        // determine minimum index of current edge
        int ind1 = _F.back()(j, i) < _F.back()((j + 1) % 3, i)
                       ? _F.back()(j, i)
                       : _F.back()((j + 1) % 3, i);
        // determine maximum index of current edge
        int ind2 = _F.back()(j, i) > _F.back()((j + 1) % 3, i)
                       ? _F.back()(j, i)
                       : _F.back()((j + 1) % 3, i);
        int l = 0;
        for (l = 0; E(ind1, l) != ind2; ++l)
          ;
        m[j] = newPind(ind1, l) + int(_P.cols());
      }
      // set new elements
      newF(0, 4 * i) = _F.back()(0, i);
      newF(1, 4 * i) = m[0];
      newF(2, 4 * i) = m[2];
      newF(0, 4 * i + 1) = m[0];
      newF(1, 4 * i + 1) = _F.back()(1, i);
      newF(2, 4 * i + 1) = m[1];
      newF(0, 4 * i + 2) = m[2];
      newF(1, 4 * i + 2) = m[1];
      newF(2, 4 * i + 2) = _F.back()(2, i);
      newF(0, 4 * i + 3) = m[0];
      newF(1, 4 * i + 3) = m[1];
      newF(2, 4 * i + 3) = m[2];
    }
    // write everything back to internal data
    _P.conservativeResize(_P.rows(), _P.cols() + newP.cols());
    _P.block(0, _P.cols() - newP.cols(), _P.rows(), newP.cols()) = newP;
    _nP.push_back(int(_P.cols()));
    _B.conservativeResize(_B.size() + newB.size());
    _B.tail(newB.size()) = newB;
    _F.push_back(newF);
    ++_maxLvl;
    return *this;
  }
  Mesh &refine(int level) {
    int currLevel = int(_F.size());
    for (auto i = currLevel; i <= level; ++i) refine();
    return *this;
  }
  void printSolution(const Eigen::VectorXd &u, const std::string &fileName) {
    std::ofstream myfile;
    auto nf = _F.back().cols();
    auto np = _PC.cols();
    myfile.open(fileName);
    // visualize the surface of the geometry in vtk format
    // vtk format header
    myfile << "# vtk DataFile Version 3.1\n";
    myfile << "this file hopefully represents my solution now\n";
    myfile << "ASCII\n";
    myfile << "DATASET UNSTRUCTURED_GRID\n";
    // print point list
    myfile << "POINTS " << _PC.cols() << " FLOAT\n";
    for (auto i = 0; i < _PC.cols(); ++i)
      myfile << std::setprecision(6) << _PC(0, i) << "\t" << _PC(1, i) << "\t"
             << u(i) << std::endl;
    // print element list
    myfile << "\nCELLS " << nf << " " << nf * 4 << std::endl;
    for (int i = 0; i < nf; ++i)
      myfile << 3 << " " << _F.back()(0, i) << " " << _F.back()(1, i) << " "
             << _F.back()(2, i) << std::endl;

    myfile << "\nCELL_TYPES " << nf << "\n";
    for (int i = 0; i < nf; ++i) myfile << 5 << std::endl;
    myfile << std::endl;

    myfile << "\nPOINT_DATA " << _PC.cols() << std::endl;
    myfile << "SCALARS Solution FLOAT\n";
    myfile << "LOOKUP_TABLE default\n";
    for (int i = 0; i < _PC.cols(); ++i) myfile << u(i) << std::endl;
    myfile << std::endl;
    myfile.close();
  }

  void printSpaceTimeGrid(const std::vector<Eigen::VectorXd> &u,
                          const Eigen::VectorXd &t,
                          const std::string &fileName) {
    std::ofstream myfile;
    auto nf = _F.back().cols();
    auto np = _PC.cols();
    myfile.open(fileName);
    // visualize the surface of the geometry in vtk format
    // vtk format header
    myfile << "# vtk DataFile Version 3.1\n";
    myfile << "this file hopefully represents my solution now\n";
    myfile << "ASCII\n";
    myfile << "DATASET UNSTRUCTURED_GRID\n";
    // print point list
    myfile << "POINTS " << _PC.cols() * t.size() << " FLOAT\n";
    for (auto j = 0; j < t.size(); ++j)
      for (auto i = 0; i < _PC.cols(); ++i)
        myfile << std::setprecision(6) << _PC(0, i) << "\t" << _PC(1, i) << "\t"
               << t(j) << std::endl;
    // print element list
    myfile << "\nCELLS " << nf * (t.size() - 1) << " "
           << nf * 7 * (t.size() - 1) << std::endl;
    for (auto j = 0; j < t.size() - 1; ++j)
      for (int i = 0; i < nf; ++i)
        myfile << 6 << " " << _F.back()(0, i) + j * np << " "
               << _F.back()(1, i) + j * np << " " << _F.back()(2, i) + j * np
               << " " << _F.back()(0, i) + (j + 1) * np << " "
               << _F.back()(1, i) + (j + 1) * np << " "
               << _F.back()(2, i) + (j + 1) * np << std::endl;

    myfile << "\nCELL_TYPES " << nf * (t.size() - 1) << "\n";
    for (int i = 0; i < nf * (t.size() - 1); ++i) myfile << 13 << std::endl;
    myfile << std::endl;

    myfile << "\nPOINT_DATA " << _PC.cols() * t.size() << std::endl;
    myfile << "SCALARS Solution FLOAT\n";
    myfile << "LOOKUP_TABLE default\n";
    for (auto j = 0; j < t.size(); ++j)
      for (int i = 0; i < _PC.cols(); ++i) myfile << u[j](i) << std::endl;
    myfile << std::endl;
    myfile.close();
  }
  void printGeometry(const Eigen::VectorXd &u, int m,
                     const std::string &fileName) {
    std::ofstream myfile;
    // double max = 0;
    int ppblock = 0;
    int **Elist = NULL;
    int M = _F.size() - 1;
    auto nfm = _F[m].cols();
    // auto nfM = _F.back().cols();
    // compute number of points per element
    ppblock = 3 * (1 << (M - m));
    Elist = (int **)calloc(1, nfm * sizeof(int *));
    for (int i = 0; i < nfm; ++i) {
      Elist[i] = (int *)calloc(1, ppblock * sizeof(int));
      // write edges
      for (int k = 0; k < 3; ++k) {
        Elist[i][k * (1 << (M - m))] = _F[m](k, i);
        if (m < M)
          recAddPoints(k, 4 * i, m + 1, &(Elist[i][k * (1 << (M - m))]));
      }
    }
    myfile.open(fileName);
    // visualize the surface of the geometry in vtk format
    // vtk format header
    myfile << "# vtk DataFile Version 3.1\n";
    myfile << "this file hopefully represents my solution now\n";
    myfile << "ASCII\n";
    myfile << "DATASET UNSTRUCTURED_GRID\n";

    // print point list
    myfile << "POINTS " << _PC.cols() << " FLOAT\n";
    for (int i = 0; i < _PC.cols(); ++i)
      myfile << std::setprecision(8) << _PC(0, i) << "\t" << _PC(1, i) << "\t"
             << u(i) << std::endl;

    // print element list
    myfile << "CELLS " << nfm << " " << nfm * (ppblock + 1) << std::endl;
    for (int i = 0; i < nfm; ++i) {
      myfile << ppblock << " ";
      for (int j = 0; j < ppblock; ++j) myfile << Elist[i][j] << " ";
      myfile << std::endl;
    }
    myfile << std::endl;

    myfile << "CELL_TYPES " << nfm << "\n";
    for (int i = 0; i < nfm; ++i) myfile << 7 << std::endl;
    myfile << std::endl;
    /* print z-values of the geometry and solved density for visualization */
    if (m == _nP.size() - 1) {
      myfile << "POINT_DATA " << _PC.cols() << std::endl;
      myfile << "SCALARS Solution FLOAT\n";
      myfile << "LOOKUP_TABLE default\n";
      for (int i = 0; i < _PC.cols(); ++i) myfile << u(i) << std::endl;
      myfile << std::endl;
    }
    for (int i = 0; i < nfm; ++i) free(Elist[i]);
    free(Elist);
    return;
  }
  unsigned int computeMaxDegree(void) const {
    auto nP = _nP[0];
    auto F = _F[0];
    Eigen::VectorXi deg(nP);
    deg.setZero();
    for (auto j = 0; j < F.cols(); ++j)
      for (auto i = 0; i < F.rows(); ++i) ++deg(F(i, j));
    return deg.maxCoeff();
  }

  void setNeumann() {
    for (auto i = 0; i < _B.size(); ++i)
      if (!_B(i)) _B(i) = 2;
  }

 private:
  /**----------------------------------------------------------------------------+
   *   int print_geometry |
   *   - visualizes the geometry and the solution in vtk format | on the
   *finest level of the discretization                                |
   *----------------------------------------------------------------------------
   **/
  void recAddPoints(int k, int i, int m, int *Elist) {
    auto M = _F.size() - 1;
    Elist[1 << (M - m)] = _F[m](k, i + 3);
    if (m < M) {
      recAddPoints(k, 4 * (i + k), m + 1, &Elist[0]);
      recAddPoints(k, 4 * (i + ((k + 1) % 3)), m + 1, &Elist[1 << (M - m)]);
    }
    return;
  }

  /**
   *  search functions for genTriang: if x = Chi.f(p,t) exists already in P
   *  then its index is returned. If not, it is attached to P and its index
   *  is returned either.
   */
  static unsigned int searchPoint(unsigned int p, const Eigen::Vector2d &t,
                                  const ParametrisationBase &Chi,
                                  Eigen::MatrixXd &P, Eigen::VectorXi &B,
                                  unsigned int *nP) {
    eigen_assert(*nP < P.cols());
    Eigen::Vector2d x = Chi.f(p, t);
    // check if point is contained in current point list
    int k = 0;
    for (k = 0; k < *nP && (P.col(k) - x).norm() > 1e-6 * P.col(k).norm(); ++k)
      ;
    // if not, attach point
    if (k == *nP) {
      P.col(*nP) = x;
      B[*nP] = Chi.boundCond(p, t);
      ++(*nP);
    }
    // return the index of the desired point
    return k;
  }
  /* initializes the coarse grid triangulation and distinguishes boundary
   * points and interior points
   */
  void genTriang(const ParametrisationBase &Chi) {
    int p = Chi.get_p();
    unsigned int nP = 0;
    Eigen::Vector2d hatp0;
    Eigen::Vector2d hatp1;
    Eigen::Vector2d hatp2;
    hatp0 << 0, 0;
    hatp1 << 1, 0;
    hatp2 << 0, 1;
    // determine worst case memory requirements
    _B.resize(3 * p);
    _B.setZero();
    _P.resize(2, 3 * p);
    _P.setZero();
    _F.clear();
    _nP.clear();
    _F.push_back(Eigen::MatrixXi());
    _F[0].resize(3, p);
    _F[0].setZero();

    // determine the vertices of the patches
    for (auto i = 0; i < p; ++i) {
      _F[0](0, i) = searchPoint(i, hatp0, Chi, _P, _B, &nP);
      _F[0](1, i) = searchPoint(i, hatp1, Chi, _P, _B, &nP);
      _F[0](2, i) = searchPoint(i, hatp2, Chi, _P, _B, &nP);
    }
    _B.conservativeResize(nP);
    _P.conservativeResize(2, nP);
    _nP.push_back(nP);
    _maxLvl = 0;
    return;
  }
  static Eigen::MatrixXi getEdgeList(const Eigen::MatrixXi &F,
                                     const Eigen::MatrixXd &P) {
    Eigen::MatrixXi E;
    E.resize(P.cols(), 7);
    E.col(0) = -Eigen::VectorXi::Ones(P.cols());
    for (auto i = 0; i < F.cols(); ++i)
      for (auto j = 0; j < 3; ++j) {
        // determine minimum index of current edge
        int ind1 = F(j, i) < F((j + 1) % 3, i) ? F(j, i) : F((j + 1) % 3, i);
        // determine maximum index of current edge
        int ind2 = F(j, i) > F((j + 1) % 3, i) ? F(j, i) : F((j + 1) % 3, i);
        // check if edge is already stored in the edge list
        int k = 0;
        for (k = 0; E(ind1, k) != -1 && E(ind1, k) != ind2; ++k)
          ;
        // if not, add it
        if (E(ind1, k) == -1) {
          // check if there is space left to store the guard -1, otherwise
          // reallocate a whole column (if this actually happens, your
          // coarse grid is crap anyway)
          if (k + 1 == E.cols()) E.conservativeResize(E.rows(), E.cols() + 1);
          E(ind1, k) = ind2;
          E(ind1, k + 1) = -1;
        }
      }
    return E;
  }
  static Eigen::MatrixXd polygonal2curved(const Eigen::MatrixXd &P,
                                          const std::vector<Eigen::MatrixXi> &F,
                                          const ParametrisationBase &Chi) {
    Eigen::MatrixXd PC = P;
    // return PC;
    Eigen::VectorXi isMapped(P.cols());
    isMapped.setZero();
    int p = int(F.front().cols());
    int n = int(F.back().cols() / F.front().cols());
    // patchwise iteration
    for (auto i = 0; i < p; ++i) {
      Eigen::Vector2d p0 = P.col(F.front()(0, i));
      Eigen::Vector2d p1 = P.col(F.front()(1, i));
      Eigen::Vector2d p2 = P.col(F.front()(2, i));
      Eigen::Vector2d a;
      for (auto j = i * n; j < (i + 1) * n; ++j)
        for (auto k = 0; k < 3; ++k)
          // check first point of the current element
          if (!isMapped(F.back()(k, j))) {
            // map point to reference element
            a = Transformations::invF_T(P.col(F.back()(k, j)), p0, p1, p2);
            // and now push it forward to actual position
            PC.col(F.back()(k, j)) = Chi.f(i, a);
            // remember that point has already been mapped
            isMapped(F.back()(k, j)) = 1;
          }
    }
    return PC;
  }
  std::vector<int> _nP;
  Eigen::VectorXi _B;
  std::vector<Eigen::MatrixXi> _F;
  Eigen::MatrixXd _P;
  Eigen::MatrixXd _PC;
  unsigned int _maxLvl;
};

}  // namespace SFEMPP
#endif
