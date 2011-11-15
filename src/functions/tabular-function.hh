#ifndef AMANZI_TABULAR_FUNCTION_HH_
#define AMANZI_TABULAR_FUNCTION_HH_

#include <vector>

#include "function.hh"

namespace Amanzi {

class TabularFunction : public Function {
 public:
  enum Form { LINEAR, CONSTANT };
 
 public:
  TabularFunction(const std::vector<double> &x, const std::vector<double> &y);
  TabularFunction(const std::vector<double> &x, const std::vector<double> &y,
                  const std::vector<Form> &form);
  ~TabularFunction() {}
  TabularFunction* Clone() const { return new TabularFunction(*this); }
  double operator() (const double *x) const;

 private:
  std::vector<double> x_, y_;
  std::vector<Form> form_;
  
 private: // helper functions
  void check_args(const std::vector<double>&, const std::vector<double>&,
     const std::vector<Form>&) const;
};

} // namespace Amanzi

#endif  // AMANZI_TABULAR_FUNCTION_HH_
