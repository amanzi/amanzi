/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (ecoon _at_ lanl.gov)
*/

//! FunctionComposition: f(x,y) = f1(x,y) * f2(x,y)

/*!

Function composition simply applies one function to the result of another.

.. math::
  f(x) = f_1( f_2(x) )

where :math:`f_1` is defined by the `"function1`" sublist, and
:math:`f_2` by the `"function2`" sublist.

* `"function1`" ``[function-spec]`` f_1 in f(x) = f_1(f_2(x))
* `"function2`" ``[function-spec]`` f_2 in f(x) = f_1(f_2(x))


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
  FunctionComposition(std::unique_ptr<Function> f1,
                      std::unique_ptr<Function> f2)
    : f1_(std::move(f1)), f2_(std::move(f2)){};
  FunctionComposition(const Function& f1, const Function& f2)
    : f1_(f1.Clone()), f2_(f2.Clone())
  {}
  FunctionComposition(const FunctionComposition& source)
    : f1_(source.f1_->Clone()), f2_(source.f2_->Clone())
  {}
  ~FunctionComposition() {} //{ if (f1_) delete f1_; if (f2_) delete f2_; }
  FunctionComposition* Clone() const { return new FunctionComposition(*this); }

  double operator()(const Kokkos::View<double*>& x) const
  {
    Kokkos::View<double*> y("y", x.extent(0));
    Kokkos::deep_copy(y, x);
    y(0) = (*f2_)(x);
    return (*f1_)(y);
  }

  void apply(const Kokkos::View<double**>& in, Kokkos::View<double*>& out) const
  {
    Kokkos::View<double*> out_1("out", out.extent(0));
    Kokkos::View<double**> tmpin("tmpin", in.extent(0), in.extent(1));
    Kokkos::deep_copy(tmpin, in);
    f2_->apply(tmpin, out_1);
    // Change all first value
    Kokkos::parallel_for(in.extent(1), KOKKOS_LAMBDA(const int& j) {
      for (int i = 0; i < in.extent(0); ++i) tmpin(i, j) = out_1(i);
    });
    f1_->apply(tmpin, out);
  }

 private:
  std::unique_ptr<Function> f1_, f2_;
  // Function *f1_, *f2_;
};

} // namespace Amanzi

#endif // AMANZI_COMPOSITION_FUNCTION_HH_
