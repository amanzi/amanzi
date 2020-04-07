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
