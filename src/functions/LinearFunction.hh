/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! LinearFunction: a multivariate linear function.

/*
  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/

/*!

A multi-variable linear function is formally defined by
 
.. math::
  f(x) = y_0 + \sum_{{j=0}}^{{n-1}} g_j (x_j - x_{{0,j}}) 

with the constant term "math:`y_0` and gradient :math:`g_0,\, g_1\,..., g_{{n-1}}`.
If the reference point :math:`x_0` is specified, it must have the same
number of values as the gradient.  Otherwise, it defaults to zero.
Note that one of the parameters in a multi-valued linear function can be time.
Here is an example:

.. code-block:: xml

  <ParameterList name="function-linear">
    <Parameter name="y0" type="double" value="1.0"/>
    <Parameter name="gradient" type="Array(double)" value="{{1.0, 2.0, 3.0}}"/>
    <Parameter name="x0" type="Array(double)" value="{{2.0, 3.0, 1.0}}"/>
  </ParameterList>

*/

#ifndef AMANZI_LINEAR_FUNCTION_HH_
#define AMANZI_LINEAR_FUNCTION_HH_

#include <vector>

#include "Function.hh"

namespace Amanzi {

class LinearFunction : public Function {
 public:
  LinearFunction(double y0, const std::vector<double> &grad);
  LinearFunction(double y0, const std::vector<double> &grad, const std::vector<double> &x0);
  ~LinearFunction() {}
  LinearFunction* Clone() const { return new LinearFunction(*this); }
  double operator()(const std::vector<double>& x) const;

 private:
  double y0_;
  std::vector<double> grad_, x0_;
};

} // namespace Amanzi

#endif
