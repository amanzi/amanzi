#ifndef AMANZI_CONSTANT_FUNCTION_HH_
#define AMANZI_CONSTANT_FUNCTION_HH_

#include "Function.hh"

namespace Amanzi {

class ConstantFunction : public Function {
 public:
  ConstantFunction(double c) : c_(c) {}
  ConstantFunction* Clone() const { return new ConstantFunction(*this); }
  double operator() (const double* x) const { return c_; }
  
 private:
  double c_;
};

} // namespace Amanzi

#endif // AMANZI_CONSTANT_FUNCTION_HH_
