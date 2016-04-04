#ifndef AMANZI_TABULAR_FUNCTION_HH_
#define AMANZI_TABULAR_FUNCTION_HH_

#include <memory>
#include <vector>

#include "Function.hh"

namespace Amanzi {

class TabularFunction : public Function {
 public:
  enum Form { LINEAR, CONSTANT, FUNCTION };
 
 public:
  TabularFunction(const std::vector<double>& x, const std::vector<double>& y,
                  const int xi);
  TabularFunction(const std::vector<double>& x, const std::vector<double>& y,
                  const int xi, const std::vector<Form>& form);
  TabularFunction(const std::vector<double>& x, const std::vector<double>& y,
                  const int xi, const std::vector<Form>& form,
                  const std::vector<Function*>& func);
  ~TabularFunction() {};

  TabularFunction* Clone() const { return new TabularFunction(*this); }
  double operator()(const std::vector<double>& x) const;

 private:
  std::vector<double> x_, y_;
  std::vector<Form> form_;
  std::vector<Function* > func_;
  int xi_;
  
 private: // helper functions
  void check_args(const std::vector<double>&, const std::vector<double>&,
     const std::vector<Form>&) const;
};

} // namespace Amanzi

#endif  // AMANZI_TABULAR_FUNCTION_HH_
