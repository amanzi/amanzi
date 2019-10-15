/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (ecoon _at_ lanl.gov)
*/

//! FunctionSeparable: f(x,y) = f1(x)*f2(y)

/*!

A separable function is defined as the product of other functions such as

.. math::
  f(x_0, x_1,...,x_{{n-1}}) = f_1(x_0)\, f_2(x_1,...,x_{{n-1}})

where :math:`f_1` is defined by the `"function1`" sublist, and
:math:`f_2` by the `"function2`" sublist.

* `"function1`" ``[function-spec]`` f_1 in f(x) = f_1(x0) * f_2(x1...)
* `"function2`" ``[function-spec]`` f_2 in f(x) = f_1(x0) * f_2(x1...)


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

class FunctionSeparable : public Function {
 public:
  FunctionSeparable(std::unique_ptr<Function> f1, std::unique_ptr<Function> f2)
    : f1_(std::move(f1)), f2_(std::move(f2)){};
  FunctionSeparable(const Function& f1, const Function& f2)
    : f1_(f1.Clone()), f2_(f2.Clone())
  {}
  FunctionSeparable(const FunctionSeparable& source)
    : f1_(source.f1_->Clone()), f2_(source.f2_->Clone())
  {}
  ~FunctionSeparable() {} //{ if (f1_) delete f1_; if (f2_) delete f2_; }
  FunctionSeparable* Clone() const { return new FunctionSeparable(*this); }

  double operator()(const Kokkos::View<double*>& x) const
  {
    // std::vector<double>::const_iterator xb = x.begin(); xb++;
    // std::vector<double> y(xb, x.end());
    auto y =
      Kokkos::subview(x,
                      Kokkos::make_pair(static_cast<size_t>(1),
                                        static_cast<size_t>(x.extent(0))));
    return (*f1_)(x) * (*f2_)(y);
  }


  void apply(const Kokkos::View<double**>& in, Kokkos::View<double*>& out) const
  {
    Kokkos::View<double**> y =
      Kokkos::subview(in,
                      Kokkos::make_pair(static_cast<size_t>(1),
                                        static_cast<size_t>(in.extent(0))),
                      Kokkos::ALL);
    Kokkos::View<double*> out_1("out_1", out.extent(0));
    f1_->apply(in, out_1);
    f2_->apply(y, out);
    Kokkos::parallel_for(out.extent(0),
                         KOKKOS_LAMBDA(const int& i) { out(i) *= out_1(i); });
  }

 private:
  std::unique_ptr<Function> f1_, f2_;
  // Function *f1_, *f2_;
};

} // namespace Amanzi

#endif // AMANZI_SEPARABLE_FUNCTION_HH_
