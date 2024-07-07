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
  : op_str_(op), amplitude_(amplitude), parameter_(parameter), shift_(shift)
{
  if (op == "cos")
    op_ = Function_kind::COS;
  else if (op == "sin")
    op_ = Function_kind::SIN;
  else if (op == "tan")
    op_ = Function_kind::TAN;
  else if (op == "acos")
    op_ = Function_kind::ACOS;
  else if (op == "asin")
    op_ = Function_kind::ASIN;
  else if (op == "atan")
    op_ = Function_kind::ATAN;
  else if (op == "cosh")
    op_ = Function_kind::COSH;
  else if (op == "sinh")
    op_ = Function_kind::SINH;
  else if (op == "tanh")
    op_ = Function_kind::TANH;
  else if (op == "exp")
    op_ = Function_kind::EXP;
  else if (op == "log")
    op_ = Function_kind::LOG;
  else if (op == "log10")
    op_ = Function_kind::LOG10;
  else if (op == "sqrt")
    op_ = Function_kind::SQRT;
  else if (op == "ceil")
    op_ = Function_kind::CEIL;
  else if (op == "fabs")
    op_ = Function_kind::FABS;
  else if (op == "abs")
    op_ = Function_kind::FABS;
  else if (op == "floor")
    op_ = Function_kind::FLOOR;
  else if (op == "mod")
    op_ = Function_kind::MOD;
  else if (op == "pow")
    op_ = Function_kind::POW;
  else if (op == "positive")
    op_ = Function_kind::POSITIVE;
  else if (op == "negative")
    op_ = Function_kind::NEGATIVE;
  else if (op == "heaviside")
    op_ = Function_kind::HEAVISIDE;
  else if (op == "sign")
    op_ = Function_kind::SIGN;
  else {
    std::stringstream m;
    m << "Invalid or unknown standard math function " << op;
    Errors::Message message(m.str());
    Exceptions::amanzi_throw(message);
  }
}


double
FunctionStandardMath::operator()(const Kokkos::View<const double**, Kokkos::HostSpace>& x) const
{
  auto f = Impl::FunctionStandardMathFunctor(op_, parameter_, amplitude_, shift_, x);
  return f(0);
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
      "FunctionStandardMath::apply1", ids_loc.extent(0), KOKKOS_LAMBDA(const int& i) {
        out(ids_loc(i)) = f(ids_loc(i));
      });
  } else {
    Kokkos::parallel_for(
      "FunctionStandardMath::apply2", in.extent(1), KOKKOS_LAMBDA(const int& i) {
        out(i) = f(i);
      });
  }
}


void
FunctionStandardMath::InvalidDomainError_(double x) const
{
  std::stringstream m;
  m << "Value " << x << " is not in the domain of operator " << op_str_ << ".";
  Errors::Message message(m.str());
  Exceptions::amanzi_throw(message);
}



} // namespace Amanzi
