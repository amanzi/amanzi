/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//! Function: base class for all functions of space and time.

/*!
Analytic, algabraic functions of space and time are used for a variety of
purposes, including boundary conditions, initial conditions, and independent
variables.

For initial conditions, functions are prescribed of space only, i.e.

:math:`u = f(x,y,z)`

For boundary conditions and independent variables, functions are also a
function of time:

:math:`u = f(t,x,y,z)`

``[function-spec]``

ONE OF:
* `"function: constant`" ``[constant-function-spec]``
OR:
* `"function: tabular`" ``[tabular-function-spec]``
OR:
* `"function: smooth step`" ``[smooth-step-function-spec]``
OR:
* `"function: polynomial`" ``[polynomial-function-spec]``
OR:
* `"function: monomial`" ``[monomial-function-spec]``
OR:
* `"function: linear`" ``[linear-function-spec]``
OR:
* `"function: separable`" ``[separable-function-spec]``
OR:
* `"function: additive`" ``[additive-function-spec]``
OR:
* `"function: multiplicative`" ``[multiplicative-function-spec]``
OR:
* `"function: composition`" ``[composition-function-spec]``
OR:
* `"function: static head`" ``[static-head-function-spec]``
OR:
* `"function: standard math`" ``[standard-math-function-spec]``
OR:
* `"function: bilinear`" ``[bilinear-function-spec]``
OR:
* `"function: distance`" ``[distance-function-spec]``
#OR:
#* `"function: squared distance`" ``[squared-distance-function-spec]``
END
*/

#ifndef AMANZI_FUNCTION_HH_
#define AMANZI_FUNCTION_HH_

#include <vector>
#include <Kokkos_Core.hpp>
#include <cassert>


namespace Amanzi {

class Function {
 public:
  virtual ~Function() {}
  virtual Function* Clone() const = 0;

  // Keep host version working
  virtual double operator()(const Kokkos::View<double*>&) const = 0;

  virtual void
  apply(const Kokkos::View<double**>&, Kokkos::View<double*>&) const = 0;
};

} // namespace Amanzi

#endif // AMANZI_FUNCTION_HH_
