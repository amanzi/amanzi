// TODO: This is very tentative.  I'm thinking it might be preferable
// to clone the constructor functions rather than hand off pointers
// to them, but cloning the polymorphic objects requires some changes
// to the base class and all the implementations.

#ifndef AMANZI_SEPARABLE_FUNCTION_HH_
#define AMANZI_SEPARABLE_FUNCTION_HH_

#include <memory>
#include <algorithm>

#include "Function.hh"

namespace Amanzi {

class SeparableFunction : public Function {
 public:
  SeparableFunction(std::unique_ptr<Function> f1, std::unique_ptr<Function> f2)
     : f1_(std::move(f1)), f2_(std::move(f2)) {};
  SeparableFunction(const Function& f1, const Function& f2)
     : f1_(f1.Clone()), f2_(f2.Clone()) {}
  SeparableFunction(const SeparableFunction& source)
     : f1_(source.f1_->Clone()), f2_(source.f2_->Clone()) {}
  ~SeparableFunction() {} //{ if (f1_) delete f1_; if (f2_) delete f2_; }
  SeparableFunction* Clone() const { return new SeparableFunction(*this); }
  double operator()(const std::vector<double>& x) const {
    std::vector<double>::const_iterator xb = x.begin(); xb++;
    std::vector<double> y(xb, x.end());
    return (*f1_)(x) * (*f2_)(y);
  }

 private:
  std::unique_ptr<Function> f1_, f2_;
  //Function *f1_, *f2_;
};

} // namespace Amanzi

#endif  // AMANZI_SEPARABLE_FUNCTION_HH_
