/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! FunctionMultiplicative: f(x,y) = f1(x,y) * f2(x,y)

/*
  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon _at_ lanl.gov)
*/

/*!

A multiplicative function simply multiplies two other function results together.

.. math::
  f(x) = f_1(x) * f_2(x)

where :math:`f_1` is defined by the `"function1`" sublist, and 
:math:`f_2` by the `"function2`" sublist.

* `"function1`" ``[function-spec]`` f_1 in f(x) = f_1(x) + f_2(x)
* `"function2`" ``[function-spec]`` f_2 in f(x) = f_1(x) + f_2(x)

Example:

.. code-block:: xml

  <ParameterList name="function-multiplicative">
    <ParameterList name="function1">
      function-specification
    </ParameterList>
    <ParameterList name="function2">
      function-specification
    </ParameterList>
  </ParameterList>
*/

#ifndef AMANZI_MULTIPLICATIVE_FUNCTION_HH_
#define AMANZI_MULTIPLICATIVE_FUNCTION_HH_

#include <memory>

#include "Function.hh"

namespace Amanzi {

class FunctionMultiplicative : public Function {
 public:
  FunctionMultiplicative(std::unique_ptr<Function> f1, std::unique_ptr<Function> f2)
     : f1_(std::move(f1)), f2_(std::move(f2)) {};
  FunctionMultiplicative(const Function& f1, const Function& f2)
     : f1_(f1.Clone()), f2_(f2.Clone()) {}
  FunctionMultiplicative(const FunctionMultiplicative& source)
     : f1_(source.f1_->Clone()), f2_(source.f2_->Clone()) {}
  ~FunctionMultiplicative() {} //{ if (f1_) delete f1_; if (f2_) delete f2_; }
  std::unique_ptr<Function> Clone() const {
    return std::make_unique<FunctionMultiplicative>(*this);
  }
  double operator()(const std::vector<double>& x) const { return (*f1_)(x) * (*f2_)(x); }

 private:
  std::unique_ptr<Function> f1_, f2_;
};

} // namespace Amanzi

#endif  // AMANZI_MULTIPLICATIVE_FUNCTION_HH_
