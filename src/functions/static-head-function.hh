#ifndef AMANZI_STATIC_HEAD_FUNCTION_HH_
#define AMANZI_STATIC_HEAD_FUNCTION_HH_

namespace Amanzi {

#include <memory>

#include "function.hh"

class StaticHeadFunction : public Function {
 public:
  StaticHeadFunction(double patm, double rho, double g, std::auto_ptr<Function> h)
      : patm_(patm), rho_g_(rho*g), h_(h) {}
  StaticHeadFunction(double patm, double rho, double g, const Function& h)
      : patm_(patm), rho_g_(rho*g), h_(h.Clone()) {}
  StaticHeadFunction(const StaticHeadFunction& src)
      : patm_(src.patm_), rho_g_(src.rho_g_), h_(src.h_->Clone()) {}
  ~StaticHeadFunction() {}
  StaticHeadFunction* Clone() const { return new StaticHeadFunction(*this); }
  // N.B. The array (t,x,y,z) is passed as *x, so that x[3] is z.
  double operator() (const double *x) const { return patm_+rho_g_*((*h_)(x)-x[3]); }
 
 private:
  double patm_, rho_g_;
  std::auto_ptr<Function> h_;
};

} // namespace Amanzi

#endif // AMANZI_STATIC_HEAD_FUNCTION_HH_
