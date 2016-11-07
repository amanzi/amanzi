/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! SeparableFunction: f(x,y) = f1(x)*f2(y)

/*
  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon _at_ lanl.gov)
*/

/*!

A separable function is defined as the product of other functions such as

.. math::
  f(x_0, x_1,...,x_{{n-1}}) = f_1(x_0)\, f_2(x_1,...,x_{{n-1}})

where :math:`f_1` is defined by the `"function1`" sublist, and 
:math:`f_2` by the `"function2`" sublist:

.. code-block:: xml

  <ParameterList name="function-separable">
    <ParameterList name="function1">
      function-specification
    </ParameterList>
    <ParameterList name="function2">
      function-specification
    </ParameterList>
  </ParameterList>

*/  

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
