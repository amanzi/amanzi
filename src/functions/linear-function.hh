#ifndef AMANZI_LINEAR_FUNCTION_HH_
#define AMANZI_LINEAR_FUNCTION_HH_

#include <vector>

#include "function.hh"

namespace Amanzi {

class LinearFunction : public Function {
 public:
  LinearFunction(double y0, const std::vector<double> &grad);
  LinearFunction(double y0, const std::vector<double> &grad, const std::vector<double> &x0);
  ~LinearFunction() {}
  LinearFunction* Clone() const { return new LinearFunction(*this); }
  double operator() (const double *x) const;

 private:
  double y0_;
  std::vector<double> grad_, x0_;
};

} // namespace Amanzi

#endif  // AMANZI_LINEAR_FUNCTION_HH_
