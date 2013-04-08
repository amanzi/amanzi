// TODO: This is very tentative.  I'm thinking it might be preferable
// to clone the constructor functions rather than hand off pointers
// to them, but cloning the polymorphic objects requires some changes
// to the base class and all the implementations.

#ifndef AMANZI_COMPOSITION_FUNCTION_HH_
#define AMANZI_COMPOSITION_FUNCTION_HH_

#include <memory>

#include "function.hh"

namespace Amanzi {

class CompositionFunction : public Function {
 public:
  CompositionFunction(std::auto_ptr<Function> f1, std::auto_ptr<Function> f2)
     : f1_(f1), f2_(f2) {};
  CompositionFunction(const Function& f1, const Function& f2)
     : f1_(f1.Clone()), f2_(f2.Clone()) {}
  CompositionFunction(const CompositionFunction& source)
     : f1_(source.f1_->Clone()), f2_(source.f2_->Clone()) {}
  ~CompositionFunction() {} //{ if (f1_) delete f1_; if (f2_) delete f2_; }
  CompositionFunction* Clone() const { return new CompositionFunction(*this); }
  double operator() (const double *x) const {
    double y = (*f2_)(x);
    return (*f1_)(&y); }

 private:
  std::auto_ptr<Function> f1_, f2_;
  //Function *f1_, *f2_;
};

} // namespace Amanzi

#endif  // AMANZI_COMPOSITION_FUNCTION_HH_
