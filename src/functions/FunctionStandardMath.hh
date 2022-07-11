/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! FunctionStandardMath: provides access to many common mathematical functions.

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

* `"operator`" ``[string]`` specifies the name of a standard mathematical
function. Available options are `"cos`", `"sin`", `"tan`", `"acos`", `"asin`",
`"atan`",
  `"cosh`", `"sinh`", `"tanh`", `"exp`", `"log`", `"log10`", `"sqrt`", `"ceil`",
  `"fabs`", `"floor`", `"mod`", and `"pow`".

* `"amplitude`" ``[double]`` specifies a multiplication factor `a` in formula `a
f(x)`. The multiplication factor is ignored by function `mod`. Default value
is 1.

* `"parameter`" ``[double]`` **1.0** specifies additional parameter `p` for
  math functions with two arguments. These functions are `"a pow(x[0], p)`"
  and `"a mod(x[0], p)`".  Alternative, scales the argument before
  application, for use in changing the period of trig functions.

* `"shift`" ``[double]`` specifies a shift of the function argument. Default is
0.

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
  enum function : int {
    COS,
    SIN,
    TAN,
    ACOS,
    ASIN,
    ATAN,
    COSH,
    SINH,
    TANH,
    EXP,
    LOG,
    LOG10,
    SQRT,
    CEIL,
    FABS,
    FLOOR,
    POW,
    MOD
  };

 public:
  FunctionStandardMath(std::string op, double amplitude, double parameter,
                       double shift);
  ~FunctionStandardMath() {}
  FunctionStandardMath* Clone() const
  {
    return new FunctionStandardMath(*this);
  }
  double operator()(const Kokkos::View<double*,Kokkos::HostSpace>&) const;

  KOKKOS_INLINE_FUNCTION double
  apply_gpu(const Kokkos::View<double**>& x, const int i) const
  {
    double x0 = x(0, i) - shift_;
    switch (op_) {
    case COS:
      return amplitude_ * cos(parameter_ * x0);
      break;
    case SIN:
      return amplitude_ * sin(parameter_ * x0);
      break;
    case TAN:
      return amplitude_ * tan(parameter_ * x0);
      break;
    case ACOS:
      return amplitude_ * acos(parameter_ * x0);
      break;
    case ASIN:
      return amplitude_ * asin(parameter_ * x0);
      break;
    case ATAN:
      return amplitude_ * atan(parameter_ * x0);
      break;
    case COSH:
      return amplitude_ * cosh(parameter_ * x0);
      break;
    case SINH:
      return amplitude_ * sinh(parameter_ * x0);
      break;
    case TANH:
      return amplitude_ * tanh(parameter_ * x0);
      break;
    case EXP:
      return amplitude_ * exp(parameter_ * x0);
      break;
    case LOG:
      assert(x0 >= 0);
      // if (x0 <= 0) InvalidDomainError_(x[0]);
      return amplitude_ * log(parameter_ * x0);
      break;
    case LOG10:
      assert(x0 >= 0);
      // if (x0 <= 0) InvalidDomainError_(x[0]);
      return amplitude_ * log10(parameter_ * x0);
      break;
    case SQRT:
      assert(x0 >= 0);
      // if (x0 < 0) InvalidDomainError_(x[0]);
      return amplitude_ * sqrt(parameter_ * x0);
      break;
    case CEIL:
      return amplitude_ * ceil(x0);
      break;
    case FABS:
      return amplitude_ * fabs(x0);
      break;
    case FLOOR:
      return amplitude_ * floor(x0);
      break;
    case POW:
      return amplitude_ * pow(x0, parameter_);
      break;
    case MOD:
      return fmod(x0, parameter_);
      break;
    default:
      printf("Invalid or unknown standard math function %d\n", op_);
      assert(false);
    }
    return 0.0;
  }

  void apply(const Kokkos::View<double**>& in, Kokkos::View<double*>& out) const
  {
    Kokkos::parallel_for(
      "FunctionStandardMath::apply",
      in.extent(1), KOKKOS_LAMBDA(const int& i) { out(i) = apply_gpu(in, i); });
  }

 private:
  void InvalidDomainError_(double x) const;

 private:
  double parameter_;
  double amplitude_;
  double shift_;
  function op_;
  // char op_[10];
};

} // namespace Amanzi

#endif // AMANZI_STANDARD_MATH_FUNCTION_HH_
