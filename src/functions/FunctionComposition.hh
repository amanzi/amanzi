/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon _at_ lanl.gov)
*/

//! FunctionComposition: f(x,y) = f1(x,y) * f2(x,y)
/*!

Function composition simply applies one function to the result of another.

.. math::
  f(x) = f_1( f_2(x) )

where :math:`f_1` is defined by the `"function1`" sublist, and
:math:`f_2` by the `"function2`" sublist.

.. _function-composition-spec:
.. admonition:: function-composition-spec

   * `"function1`" ``[function-spec]`` :math:`f_1` in :math:`f(x) = f_1(f_2(x))`
   * `"function2`" ``[function-spec]`` :math:`f_2` in :math:`f(x) = f_1(f_2(x))`


.. code-block:: xml

  <ParameterList name="function-composition">
    <ParameterList name="function1">
      function-specification
    </ParameterList>
    <ParameterList name="function2">
      function-specification
    </ParameterList>
  </ParameterList>
*/

#ifndef AMANZI_COMPOSITION_FUNCTION_HH_
#define AMANZI_COMPOSITION_FUNCTION_HH_

#include <memory>

#include "Function.hh"

namespace Amanzi {

class FunctionComposition : public Function {
 public:
  FunctionComposition(std::unique_ptr<Function> f1, std::unique_ptr<Function> f2)
    : f1_(std::move(f1)), f2_(std::move(f2)){};
  FunctionComposition(const Function& f1, const Function& f2) : f1_(f1.Clone()), f2_(f2.Clone()) {}
  FunctionComposition(const FunctionComposition& source)
    : f1_(source.f1_->Clone()), f2_(source.f2_->Clone())
  {}
  ~FunctionComposition() {} //{ if (f1_) delete f1_; if (f2_) delete f2_; }
  std::unique_ptr<Function> Clone() const { return std::make_unique<FunctionComposition>(*this); }
  double operator()(const std::vector<double>& x) const
  {
    std::vector<double> y(x);
    y[0] = (*f2_)(x);
    return (*f1_)(y);
  }

 private:
  std::unique_ptr<Function> f1_, f2_;
};

} // namespace Amanzi

#endif // AMANZI_COMPOSITION_FUNCTION_HH_
