#ifndef AMANZI_POINTER_FUNCTION_HH_
#define AMANZI_POINTER_FUNCTION_HH_

#include <vector>

#include "function.hh"

namespace Amanzi {

class PointerFunction : public Function {
 public:
  PointerFunction(double (*f)(const double*, const double*)) : f_(f), np_(0), p_(0) {}
  PointerFunction(double (*f)(const double*, const double*), const std::vector<double>&);
  PointerFunction(const PointerFunction&);
  ~PointerFunction();
  PointerFunction* Clone() const { return new PointerFunction(*this); }
  double operator() (const double *x) const { return (*f_)(x, p_); }
  
 private:
  double (*f_)(const double*, const double*);
  int np_;
  double *p_;
};

} // namespace Amanzi

#endif // AMANZI_POINTER_FUNCTION_HH_
