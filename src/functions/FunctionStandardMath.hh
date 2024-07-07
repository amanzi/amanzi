/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
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

#pragma once

#include <vector>
#include <string>

#include "Function.hh"

namespace Amanzi {

enum class Function_kind : int {
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
  MOD,
  POSITIVE,
  NEGATIVE,
  HEAVISIDE,
  SIGN
};

class FunctionStandardMath : public Function {
 public:
  FunctionStandardMath(std::string op, double amplitude, double parameter, double shift);

  std::unique_ptr<Function> Clone() const override { return std::make_unique<FunctionStandardMath>(*this); }

  double operator()(const Kokkos::View<const double**, Kokkos::HostSpace>&) const override;

  void apply(const Kokkos::View<const double**>& in,
             Kokkos::View<double*>& out,
             const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const override;

 private:
  void InvalidDomainError_(double x) const;

 private:
  std::string op_str_;
  Function_kind op_;
  double parameter_;
  double amplitude_;
  double shift_;
};



namespace Impl {

template <class InView_type>
class FunctionStandardMathFunctor {
 public:
  FunctionStandardMathFunctor(const Function_kind& op,
          double parameter,
          double amplitude,
          double shift,
          const InView_type& in)
    : op_(op), parameter_(parameter), amplitude_(amplitude), shift_(shift), in_(in) {}

  KOKKOS_INLINE_FUNCTION
  double operator()(const int i) const
  {
    double x0 = in_(0, i) - shift_;
    switch (op_) {
    case Function_kind::COS:
      return amplitude_ * Kokkos::cos(parameter_ * x0);
      break;
    case Function_kind::SIN:
      return amplitude_ * Kokkos::sin(parameter_ * x0);
      break;
    case Function_kind::TAN:
      return amplitude_ * Kokkos::tan(parameter_ * x0);
      break;
    case Function_kind::ACOS:
      return amplitude_ * Kokkos::acos(parameter_ * x0);
      break;
    case Function_kind::ASIN:
      return amplitude_ * Kokkos::asin(parameter_ * x0);
      break;
    case Function_kind::ATAN:
      return amplitude_ * Kokkos::atan(parameter_ * x0);
      break;
    case Function_kind::COSH:
      return amplitude_ * Kokkos::cosh(parameter_ * x0);
      break;
    case Function_kind::SINH:
      return amplitude_ * Kokkos::sinh(parameter_ * x0);
      break;
    case Function_kind::TANH:
      return amplitude_ * Kokkos::tanh(parameter_ * x0);
      break;
    case Function_kind::EXP:
      return amplitude_ * Kokkos::exp(parameter_ * x0);
      break;
    case Function_kind::LOG:
      assert(x0 >= 0);
      // if (x0 <= 0) InvalidDomainError_(x[0]);
      return amplitude_ * Kokkos::log(parameter_ * x0);
      break;
    case Function_kind::LOG10:
      assert(x0 >= 0);
      // if (x0 <= 0) InvalidDomainError_(x[0]);
      return amplitude_ * Kokkos::log10(parameter_ * x0);
      break;
    case Function_kind::SQRT:
      assert(x0 >= 0);
      // if (x0 < 0) InvalidDomainError_(x[0]);
      return amplitude_ * Kokkos::sqrt(parameter_ * x0);
      break;
    case Function_kind::CEIL:
      return amplitude_ * Kokkos::ceil(x0);
      break;
    case Function_kind::FABS:
      return amplitude_ * Kokkos::abs(x0);
      break;
    case Function_kind::FLOOR:
      return amplitude_ * Kokkos::floor(x0);
      break;
    case Function_kind::POW:
      return amplitude_ * Kokkos::pow(x0, parameter_);
      break;
    case Function_kind::MOD:
      return Kokkos::fmod(x0, parameter_);
      break;
    case Function_kind::POSITIVE:
      return amplitude_ * (x0 > 0 ? x0 : 0);
      break;
    case Function_kind::NEGATIVE:
      return amplitude_ * (x0 < 0 ? x0 : 0);
      break;
    case Function_kind::HEAVISIDE:
      return amplitude_ * (x0 > 0 ? 1 : 0);
      break;
    case Function_kind::SIGN:
      return amplitude_ * (x0 > 0 ? 1 : (x0 < 0 ? -1 : 0));
      break;
    default:
      assert(false);
    }
    return 0.0;
  }

 private:
  Function_kind op_;
  double parameter_;
  double amplitude_;
  double shift_;

  InView_type in_;
};

} // namespace Impl
} // namespace Amanzi

