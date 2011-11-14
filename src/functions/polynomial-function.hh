#ifndef AMANZI_POLYNOMIAL_FUNCTION_HH_
#define AMANZI_POLYNOMIAL_FUNCTION_HH_

#include <vector>

#include "function.hh"

namespace Amanzi {

class PolynomialFunction : public Function {
 public:
  PolynomialFunction(const std::vector<double> &c, const std::vector<int> &p, double x0 = 0.0);
  ~PolynomialFunction() {}
  PolynomialFunction* Clone() const { return new PolynomialFunction(*this); }
  double operator() (const double *x) const;

 private:
  int pmin_;
  int pmax_;
  double x0_;
  std::vector<double> c_;
};

} // namespace Amanzi

#endif  // AMANZI_POLYNOMIAL_FUNCTION_HH_
