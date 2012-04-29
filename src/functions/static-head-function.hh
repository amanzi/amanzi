#ifndef AMANZI_STATIC_HEAD_FUNCTION_HH_
#define AMANZI_STATIC_HEAD_FUNCTION_HH_

namespace Amanzi {

#include <memory>

#include "function.hh"

class StaticHeadFunction : public Function {
 public:
  StaticHeadFunction(double patm, double rho, double g, std::auto_ptr<Function> h, int dim)
      : patm_(patm), rho_g_(rho*g), h_(h), dim_(dim) {}
  StaticHeadFunction(double patm, double rho, double g, const Function& h, int dim)
      : patm_(patm), rho_g_(rho*g), h_(h.Clone()), dim_(dim) {}
  StaticHeadFunction(const StaticHeadFunction& src)
      : patm_(src.patm_), rho_g_(src.rho_g_), h_(src.h_->Clone()), dim_(src.dim_) {}
  ~StaticHeadFunction() {}
  StaticHeadFunction* Clone() const { return new StaticHeadFunction(*this); }
  // The array (t,x,y,z) is passed as *x, so that x[dim_] is z in 3D, y in 2D
  double operator() (const double *x) const { return patm_+rho_g_ * ((*h_)(x)-x[dim_]); }
 
 private:
  int dim_;
  double patm_, rho_g_;
  std::auto_ptr<Function> h_;
};

} // namespace Amanzi

#endif // AMANZI_STATIC_HEAD_FUNCTION_HH_
