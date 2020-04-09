/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include <math.h>

#include "errors.hh"
#include "FunctionMonomial.hh"

namespace Amanzi {

FunctionMonomial::FunctionMonomial(double c, const Kokkos::View<double*>& x0,
                                   const Kokkos::View<int*>& p)
{
  if (x0.extent(0) != p.extent(0)) {
    Errors::Message m;
    m << "Mismatch of multi-index and reference point dimensions.";
    Exceptions::amanzi_throw(m);
  }
  c_ = c;
  x0_ = x0;
  p_ = p;
}

double
FunctionMonomial::operator()(const Kokkos::View<double*>& x) const
{
  double y = c_;
  if (x.extent(0) < x0_.extent(0)) {
    Errors::Message m;
    m << "FunctionMonomial expects higher-dimensional argument.";
    Exceptions::amanzi_throw(m);
  }
  for (int j = 0; j < x0_.extent(0); ++j) y *= pow(x[j] - x0_[j], p_[j]);
  return y;
}

} // namespace Amanzi
