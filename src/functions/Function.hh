/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! Function: base class for all functions of space and time.

/*
  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/

/*!
Analytic, algabraic functions of space and time are used for a variety of
purposes, including boundary conditions, initial conditions, and independent
variables.

For initial conditions, functions are prescribed of space only, i.e.

:math:`u = f(x,y,z)`

For boundary conditions and independent variables, functions are also a
function of time:

:math:`u = f(t,x,y,z)`

Note, this does not follow the `"typed`" format for legacy reasons.

.. function-spec:
.. admonition:: function-spec

  ONE OF:

  * `"function: constant`" ``[function-constant-spec]``

  OR:

  * `"function: tabular`" ``[function-tabular-spec]``

  OR:

  * `"function: smooth step`" ``[function-smooth-step-spec]``

  OR:

  * `"function: polynomial`" ``[function-polynomial-spec]``

  OR:

  * `"function: monomial`" ``[function-monomial-spec]``

  OR:

  * `"function: linear`" ``[function-linear-spec]``

  OR:

  * `"function: separable`" ``[function-separable-spec]``

  OR:

  * `"function: additive`" ``[function-additive-spec]``

  OR:

  * `"function: multiplicative`" ``[function-multiplicative-spec]``

  OR:

  * `"function: composition`" ``[function-composition-spec]``

  OR:

  * `"function: static head`" ``[function-static-head-spec]``

  OR:

  * `"function: standard math`" ``[function-standard-math-spec]``

  OR:

  * `"function: bilinear`" ``[function-bilinear-spec]``

  OR:

  * `"function: distance`" ``[function-distance-spec]``

  END

*/
  
#ifndef AMANZI_FUNCTION_HH_
#define AMANZI_FUNCTION_HH_

#include <vector>

namespace Amanzi {

class Function {
 public:
  virtual ~Function() {}
  virtual Function* Clone() const = 0;
  virtual double operator()(const std::vector<double>& ) const = 0;
};

} // namespace Amanzi

#endif // AMANZI_FUNCTION_HH_
