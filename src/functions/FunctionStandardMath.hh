/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Provides access to many common mathematical functions.
/*!

These functions allow to set up non-trivial time-dependent boundary conditions
which increases a set of analytic solutions that can be used in convergence
analysis tests.

.. math::
  f(x) = A * operator( p * (x - s) )

or

.. math::
  f(x) = A * operator(x-s, p)

Note that these operate only on the first coordinate, which is often time.
Function composition can be used to apply these to other coordinates (or
better yet a dimension could/should be added upon request).

.. _function-standard-math-spec:
.. admonition:: function-standard-math-spec

   * `"operator`" ``[string]`` specifies the name of a standard mathematical
     function.  Available options are:

     - trigonometric operators: `"cos`", `"sin`", `"tan`", `"acos`", `"asin`",
       `"atan`"
     - hyperbolic trig operators: `"cosh`", `"sinh`", `"tanh`"
     - power/log operators: `"pow`", `"exp`", `"log`", `"log10`", `"sqrt`",
     - integral operators: `"ceil`", `"floor`", `"mod`",
     - `"abs`", `"fabs`", `"positive`" (0 for negative values), `"negative`" (0
       for positive values), `"heaviside`", `"sign`"

   * `"amplitude`" ``[double]`` specifies a multiplication factor `a` in
     formula `a f(x)`.  The multiplication factor is ignored by function
     `mod`. Default value is 1.

   * `"parameter`" ``[double]`` **1.0** specifies additional parameter `p` for
     math functions with two arguments. These functions are `"a pow(x[0], p)`"
     and `"a mod(x[0], p)`".  Alternative, scales the argument before
     application, for use in changing the period of trig functions.

   * `"shift`" ``[double]`` specifies a shift of the function argument. Default
     is 0.

Example:

.. code-block:: xml

  <ParameterList name="function-standard-math">
    <Parameter name="operator" type="string" value="sqrt"/>
    <Parameter name="amplitude" type="double" value="1e-7"/>
    <Parameter name="shift" type="double" value="0.1"/>
  </ParameterList>

This example defines function `1e-7 sqrt(t-0.1)`.
 */

#ifndef AMANZI_STANDARD_MATH_FUNCTION_HH_
#define AMANZI_STANDARD_MATH_FUNCTION_HH_

#include <vector>
#include <string>

#include "Function.hh"

namespace Amanzi {

class FunctionStandardMath : public Function {
 public:
  FunctionStandardMath(std::string op, double amplitude, double parameter, double shift);
  ~FunctionStandardMath() {}
  std::unique_ptr<Function> Clone() const { return std::make_unique<FunctionStandardMath>(*this); }
  double operator()(const std::vector<double>& x) const;

 private:
  void InvalidDomainError_(double x) const;

 private:
  double parameter_;
  double amplitude_;
  double shift_;
  std::string op_;
};

} // namespace Amanzi

#endif // AMANZI_STANDARD_MATH_FUNCTION_HH_
