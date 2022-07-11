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

FunctionMonomial::FunctionMonomial(double c, const Kokkos::View<double*,Kokkos::HostSpace>& x0,
                                   const Kokkos::View<int*,Kokkos::HostSpace>& p)
{
  if (x0.extent(0) != p.extent(0)) {
    Errors::Message m;
    m << "Mismatch of multi-index and reference point dimensions.";
    Exceptions::amanzi_throw(m);
  }
  c_ = c;
  Kokkos::resize(x0_,x0.extent(0)); 
  Kokkos::resize(p_,p.extent(0));
  Kokkos::deep_copy(x0_.view_host(),x0);
  Kokkos::deep_copy(p_.view_host(),p);
  Kokkos::deep_copy(x0_.view_device(),x0);
  Kokkos::deep_copy(p_.view_device(),p);   
}

double
FunctionMonomial::operator()(const Kokkos::View<double*,Kokkos::HostSpace>& x) const
{
  double y = c_;
  if (x.extent(0) < x0_.extent(0)) {
    Errors::Message m;
    m << "FunctionMonomial expects higher-dimensional argument.";
    Exceptions::amanzi_throw(m);
  }
  for (int j = 0; j < x0_.extent(0); ++j) y *= pow(x[j] - x0_.view_host()[j], p_.view_host()[j]);
  return y;
}

} // namespace Amanzi
