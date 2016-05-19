/*
  Functions

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  TODO: This is very tentative.  I'm thinking it might be preferable
  to clone the constructor functions rather than hand off pointers
  to them, but cloning the polymorphic objects requires some changes
 to the base class and all the implementations.
*/

#ifndef AMANZI_ADDITIVE_FUNCTION_HH_
#define AMANZI_ADDITIVE_FUNCTION_HH_

#include <memory>

#include "Function.hh"

namespace Amanzi {

class AdditiveFunction : public Function {
 public:
  AdditiveFunction(std::auto_ptr<Function> f1, std::auto_ptr<Function> f2)
     : f1_(f1), f2_(f2) {};
  AdditiveFunction(const Function& f1, const Function& f2)
     : f1_(f1.Clone()), f2_(f2.Clone()) {}
  AdditiveFunction(const AdditiveFunction& source)
     : f1_(source.f1_->Clone()), f2_(source.f2_->Clone()) {}
  ~AdditiveFunction() {} //{ if (f1_) delete f1_; if (f2_) delete f2_; }
  AdditiveFunction* Clone() const { return new AdditiveFunction(*this); }
  double operator()(const std::vector<double>& x) const { return (*f1_)(x) + (*f2_)(x); }

 private:
  std::auto_ptr<Function> f1_, f2_;
  //Function *f1_, *f2_;
};

} // namespace Amanzi

#endif 
