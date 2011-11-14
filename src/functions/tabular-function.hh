#ifndef AMANZI_TABULAR_FUNCTION_HH_
#define AMANZI_TABULAR_FUNCTION_HH_

#include <vector>

#include "function.hh"

namespace Amanzi {

class TabularFunction : public Function {
 public:
  TabularFunction(const std::vector<double> &x, const std::vector<double> &y);
  ~TabularFunction() {}
  TabularFunction* Clone() const { return new TabularFunction(*this); }
  double operator() (const double *x) const;

 private:
  std::vector<double> x_, y_;
};

} // namespace Amanzi

#endif  // AMANZI_TABULAR_FUNCTION_HH_
