//
// Copyright (c) 2023, Michael Multerer
//
// All rights reserved.
//
// This source code is subject to the GNU Affero General Public License v3.0
// license and without any warranty, see <https://github.com/muchip/FMCA>
// for further information.
//

#ifndef GEOMETRY_PARAMETRISATION_H_
#define GEOMETRY_PARAMETRISATION_H_

namespace SFEMPP {
/**
 *  patch structure
 **/
struct Patch {
  constexpr Patch(Eigen::Vector2d (*_f)(const Eigen::Vector2d &),
                  Eigen::Vector2d (*_dx)(const Eigen::Vector2d &),
                  Eigen::Vector2d (*_dy)(const Eigen::Vector2d &),
                  double (*_det)(const Eigen::Vector2d &))
      : f(_f), dx(_dx), dy(_dy), det(_det) {}
  Eigen::Vector2d (*f)(const Eigen::Vector2d &);
  Eigen::Vector2d (*dx)(const Eigen::Vector2d &);
  Eigen::Vector2d (*dy)(const Eigen::Vector2d &);
  double (*det)(const Eigen::Vector2d &);
};

/**
 *  patch members: all of them have (0,0) as a vertex and a curve parametrised
 *  over [0,1]. Curve has to provide a .f and a .df method
 **/
template <unsigned int P, typename Curve>
Eigen::Vector2d Chi(const Eigen::Vector2d &a) {
  return Transformations::tildeF_T(a, Eigen::Vector2d::Zero(), Curve(P));
}
template <unsigned int P, typename Curve>
Eigen::Vector2d dChidx(const Eigen::Vector2d &a) {
  return Transformations::dtildeF_Tdx(a, Eigen::Vector2d::Zero(), Curve(P));
}
template <unsigned int P, typename Curve>
Eigen::Vector2d dChidy(const Eigen::Vector2d &a) {
  return Transformations::dtildeF_Tdy(a, Eigen::Vector2d::Zero(), Curve(P));
}
template <unsigned int P, typename Curve>
double det_dChi(const Eigen::Vector2d &a) {
  return Transformations::det_dtildeF_T(a, Eigen::Vector2d::Zero(), Curve(P));
}

/**
 *   a compile time initialised array to wrap the patches
 **/
template <unsigned int P, typename T>
struct PatchVector {
  static constexpr Patch currPatch = Patch(
      Chi<P - 1, T>, dChidx<P - 1, T>, dChidy<P - 1, T>, det_dChi<P - 1, T>);
  PatchVector<P - 1, T> otherPatches;
  constexpr const Patch &operator[](int const &i) const;
};

template <typename T>
struct PatchVector<1, T> {
  static constexpr Patch currPatch =
      Patch(Chi<0, T>, dChidx<0, T>, dChidy<0, T>, det_dChi<0, T>);
  constexpr const Patch &operator[](int const &i) const;
};

template <unsigned int P, typename Curve>
class Parametrisation : public ParametrisationBase {
 public:
  unsigned int boundCond(unsigned int p, const Eigen::Vector2d &a) const {
    (void)p;
    if (abs(a(0) + a(1) - 1) < SFEMPP_ZERO_TOLERANCE)
      return (0);
    else
      return 1;
  }
  unsigned int get_p() const { return p_; }
  Eigen::Vector2d f(unsigned int p, const Eigen::Vector2d &a) const {
    return Patches_[p].f(a);
  }
  Eigen::Vector2d dx(unsigned int p, const Eigen::Vector2d &a) const {
    return Patches_[p].dx(a);
  }
  Eigen::Vector2d dy(unsigned int p, const Eigen::Vector2d &a) const {
    return Patches_[p].dy(a);
  }
  double det(unsigned int p, const Eigen::Vector2d &a) const {
    return Patches_[p].det(a);
  }

 protected:
  unsigned p_ = P;
  PatchVector<P, Curve> Patches_;
};

// necessary defininitions
template <unsigned int P, typename T>
constexpr const Patch &PatchVector<P, T>::operator[](int const &i) const {
  return (i == (P - 1) ? currPatch : otherPatches[i]);
}
template <typename T>
constexpr const Patch &PatchVector<1, T>::operator[](int const &i) const {
  return currPatch;
}
template <unsigned int P, typename T>
constexpr Patch PatchVector<P, T>::currPatch;
template <typename T>
constexpr Patch PatchVector<1, T>::currPatch;
}  // namespace SFEMPP
#endif
