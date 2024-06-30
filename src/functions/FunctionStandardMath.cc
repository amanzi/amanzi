/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//!
#include <math.h>

#include "FunctionStandardMath.hh"
#include "errors.hh"

namespace Amanzi {

FunctionStandardMath::FunctionStandardMath(std::string op,
                                           double amplitude,
                                           double parameter,
                                           double shift)
  : amplitude_(amplitude), parameter_(parameter), shift_(shift)
{
  if (op == "cos")
    op_ = COS;
  else if (op == "sin")
    op_ = SIN;
  else if (op == "tan")
    op_ = TAN;
  else if (op == "acos")
    op_ = ACOS;
  else if (op == "asin")
    op_ = ASIN;
  else if (op == "atan")
    op_ = ATAN;
  else if (op == "cosh")
    op_ = COSH;
  else if (op == "sinh")
    op_ = SINH;
  else if (op == "tanh")
    op_ = TANH;
  else if (op == "exp")
    op_ = EXP;
  else if (op == "log")
    op_ = LOG;
  else if (op == "log10")
    op_ = LOG10;
  else if (op == "sqrt")
    op_ = SQRT;
  else if (op == "ceil")
    op_ = CEIL;
  else if (op == "fabs")
    op_ = FABS;
  else if (op == "abs")
    op_ = FABS;
  else if (op == "floor")
    op_ = FLOOR;
  else if (op == "mod")
    op_ = MOD;
  else if (op == "pow")
    op_ = POW;
  else if (op == "positive")
    op_ = POSITIVE;
  else if (op == "negative")
    op_ = NEGATIVE;
  else if (op == "heaviside")
    op_ = HEAVISIDE;
  else if (op == "sign")
    op_ = SIGN;
  else {
    std::stringstream m;
    m << "Invalid or unknown standard math function " << op;
    Errors::Message message(m.str());
    Exceptions::amanzi_throw(message);
  }
}


double
FunctionStandardMath::operator()(const Kokkos::View<double*, Kokkos::HostSpace>& x) const
{
  double x0 = x[0] - shift_;
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
    if (x0 <= 0) InvalidDomainError_(x[0]);
    return amplitude_ * log(parameter_ * x0);
    break;
  case LOG10:
    if (x0 <= 0) InvalidDomainError_(x[0]);
    return amplitude_ * log10(parameter_ * x0);
    break;
  case SQRT:
    if (x0 < 0) InvalidDomainError_(x[0]);
    return amplitude_ * sqrt(parameter_ * x0);
    break;
  case CEIL:
    return amplitude_ * ceil(x0);
    break;
  case FABS:
    return amplitude_ * fabs(x0);
  case FLOOR:
    return amplitude_ * floor(x0);
    break;
  case POW:
    return amplitude_ * pow(x0, parameter_);
    break;
  case MOD:
    return fmod(x0, parameter_);
    break;
  case POSITIVE:
    return amplitude_ * (x0 > 0 ? x0 : 0);
    break;
  case NEGATIVE:
    return amplitude_ * (x0 < 0 ? x0 : 0);
    break;
  case HEAVISIDE:
    return amplitude_ * (x0 > 0 ? 1 : 0);
    break;
  case SIGN:
    return amplitude_ * (x0 > 0 ? 1 : (x0 < 0 ? -1 : 0));
    break;
  default:
    std::stringstream m;
    m << "Invalid or unknown standard math function " << op_;
    Errors::Message message(m.str());
    Exceptions::amanzi_throw(message);
  }
  return 0.0;
}


void
FunctionStandardMath::apply(const Kokkos::View<const double**>& in,
                            Kokkos::View<double*>& out,
                            const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const
{
  auto f = Impl::FunctionStandardMathFunctor(op_, parameter_, amplitude_, shift_, in);
  if (ids) {
    auto ids_loc = *ids;
    Kokkos::parallel_for(
      "FunctionStandardMath::apply1", ids_loc.extent(0), KOKKOS_CLASS_LAMBDA(const int& i) {
        out(ids_loc(i)) = f(ids_loc(i));
      });
  } else {
    Kokkos::parallel_for(
      "FunctionStandardMath::apply2", in.extent(1), KOKKOS_CLASS_LAMBDA(const int& i) {
        out(i) = f(i);
      });
  }
}


void
FunctionStandardMath::InvalidDomainError_(double x) const
{
  std::stringstream m;
  m << "Value " << x << " is not in the domain of operator " << op_ << ".";
  Errors::Message message(m.str());
  Exceptions::amanzi_throw(message);
}



} // namespace Amanzi
