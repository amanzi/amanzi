#ifndef AMANZI_SMOOTH_STEP_FUNCTION_HH_
#define AMANZI_SMOOTH_STEP_FUNCTION_HH_

#include "function.hh"

namespace Amanzi {

class SmoothStepFunction : public Function {
 public:
  SmoothStepFunction(double x0, double y0, double x1, double y1);
  ~SmoothStepFunction() {};
  SmoothStepFunction* Clone() const { return new SmoothStepFunction(*this); }
  double operator() (const double *x) const;

 private:
  double x0_, y0_, x1_, y1_;
};

} // namespace Amanzi

#endif  // AMANZI_SMOOTH_STEP_FUNCTION_HH_
