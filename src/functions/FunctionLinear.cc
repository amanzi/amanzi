/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "FunctionLinear.hh"
#include "errors.hh"

namespace Amanzi {

FunctionLinear::FunctionLinear(double y0, const Kokkos::View<double*,Kokkos::HostSpace>& grad,
                               const Kokkos::View<double*,Kokkos::HostSpace>& x0)
{
  if (grad.extent(0) < 1) {
    Errors::Message m;
    m << "At least one value required for the gradient vector";
    Exceptions::amanzi_throw(m);
  }
  if (x0.extent(0) != grad.extent(0)) {
    Errors::Message m;
    m << "Mismatch of gradient and point dimensions.";
    Exceptions::amanzi_throw(m);
  }
  y0_ = y0;

  Kokkos::resize(x0_,x0.extent(0)); 
  Kokkos::resize(grad_,grad.extent(0)); 

  Kokkos::deep_copy(x0_.view_host(),x0); 
  Kokkos::deep_copy(grad_.view_host(),grad); 

  Kokkos::deep_copy(x0_.view_device(),x0); 
  Kokkos::deep_copy(grad_.view_device(),grad); 
}

FunctionLinear::FunctionLinear(double y0, const Kokkos::View<double*,Kokkos::HostSpace>& grad)
{
  if (grad.extent(0) < 1) {
    Errors::Message m;
    m << "at least one value required for the gradient vector";
    Exceptions::amanzi_throw(m);
  }
  y0_ = y0;

  Kokkos::resize(x0_, grad.extent(0));
  for (int i = 0; i < x0_.extent(0); ++i) { x0_.view_host()(i) = 0.0; }

  Kokkos::resize(grad_,grad.extent(0)); 
  Kokkos::deep_copy(grad_.view_host(),grad); 
  Kokkos::deep_copy(x0_.view_device(),x0_.view_host()); 
  Kokkos::deep_copy(grad_.view_device(),grad); 
}

double
FunctionLinear::operator()(const Kokkos::View<double*,Kokkos::HostSpace>& x) const
{
  double y = y0_;
  if (x.extent(0) < grad_.extent(0)) {
    Errors::Message m;
    m << "FunctionLinear expects higher-dimensional argument.";
    Exceptions::amanzi_throw(m);
  }
  for (int j = 0; j < grad_.extent(0); ++j) y += grad_.view_host()[j] * (x[j] - x0_.view_host()[j]);
  return y;
}

} // namespace Amanzi
