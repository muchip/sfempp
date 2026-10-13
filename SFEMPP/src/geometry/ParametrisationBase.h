#ifndef __PARAMETRISATION__CLASS__
#define __PARAMETRISATION__CLASS__

namespace SFEMPP {
class ParametrisationBase {
 public:
  virtual unsigned int boundCond(unsigned int p,
                                 const Eigen::Vector2d &a) const = 0;
  virtual unsigned int get_p() const = 0;
  virtual Eigen::Vector2d f(unsigned int p, const Eigen::Vector2d &a) const = 0;
  virtual Eigen::Vector2d dx(unsigned int p,
                             const Eigen::Vector2d &a) const = 0;
  virtual Eigen::Vector2d dy(unsigned int p,
                             const Eigen::Vector2d &a) const = 0;
  virtual double det(unsigned int p, const Eigen::Vector2d &a) const = 0;
};
}  // namespace SFEMPP
#endif
