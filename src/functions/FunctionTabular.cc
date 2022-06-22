/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "FunctionTabular.hh"
#include "errors.hh"

namespace Amanzi {

FunctionTabular::FunctionTabular(const Kokkos::View<double*,Kokkos::HostSpace>& x,
                                 const Kokkos::View<double*,Kokkos::HostSpace>& y, const int xi)
  : xi_(xi)
{
  Kokkos::resize(x_,x.extent(0)); 
  Kokkos::resize(y_,y.extent(0)); 

  Kokkos::deep_copy(x_.view_host(),x);
  Kokkos::deep_copy(y_.view_host(),y);

  Kokkos::resize(form_, x.extent(0) - 1);
  for (int i = 0; i < form_.extent(0); ++i) { form_.view_host()(i) = LINEAR; }
  check_args(x, y, form_.view_host());

  Kokkos::deep_copy(x_.view_device(),x_.view_host());
  Kokkos::deep_copy(y_.view_device(),y_.view_host());
  Kokkos::deep_copy(form_.view_device(),form_.view_host());

}

FunctionTabular::FunctionTabular(const Kokkos::View<double*,Kokkos::HostSpace>& x,
                                 const Kokkos::View<double*,Kokkos::HostSpace>& y, const int xi,
                                 const Kokkos::View<Form*,Kokkos::HostSpace>& form)
  : xi_(xi)
{
  Kokkos::resize(x_,x.extent(0)); 
  Kokkos::resize(y_,y.extent(0)); 
  Kokkos::resize(form_,form.extent(0)); 

  Kokkos::deep_copy(x_.view_host(),x);
  Kokkos::deep_copy(y_.view_host(),y);
  Kokkos::deep_copy(form_.view_host(),form);
  check_args(x, y, form);


  Kokkos::deep_copy(x_.view_device(),x_.view_host());
  Kokkos::deep_copy(y_.view_device(),y_.view_host());
  Kokkos::deep_copy(form_.view_device(),form_.view_host());

}

FunctionTabular::FunctionTabular(const Kokkos::View<double*,Kokkos::HostSpace>& x,
                                 const Kokkos::View<double*,Kokkos::HostSpace>& y, const int xi,
                                 const Kokkos::View<Form*,Kokkos::HostSpace>& form,
                                 const std::vector<Function*>& func)
  : xi_(xi), func_(func)
{
  Kokkos::resize(x_,x.extent(0)); 
  Kokkos::resize(y_,y.extent(0)); 
  Kokkos::resize(form_,form.extent(0)); 
  Kokkos::deep_copy(x_.view_host(),x);
  Kokkos::deep_copy(y_.view_host(),y);
  Kokkos::deep_copy(form_.view_host(),form);
  check_args(x, y, form);
  Kokkos::deep_copy(x_.view_device(),x_.view_host());
  Kokkos::deep_copy(y_.view_device(),y_.view_host());
  Kokkos::deep_copy(form_.view_device(),form_.view_host());
}

void
FunctionTabular::check_args(const Kokkos::View<double*,Kokkos::HostSpace>& x,
                            const Kokkos::View<double*,Kokkos::HostSpace>& y,
                            const Kokkos::View<Form*,Kokkos::HostSpace>& form) const
{
  if (x.extent(0) != y.extent(0)) {
    Errors::Message m;
    m << "the number of x and y values differ";
    Exceptions::amanzi_throw(m);
  }
  if (x.extent(0) < 2) {
    Errors::Message m;
    m << "at least two table values must be given";
    Exceptions::amanzi_throw(m);
  }
  for (int j = 1; j < x.extent(0); ++j) {
    if (x[j] <= x[j - 1]) {
      Errors::Message m;
      m << "x values are not strictly increasing";
      Exceptions::amanzi_throw(m);
    }
  }
  if (form.extent(0) != x.extent(0) - 1) {
    Errors::Message m;
    m << "incorrect number of form values specified";
    Exceptions::amanzi_throw(m);
  }
}

double
FunctionTabular::operator()(const Kokkos::View<double*,Kokkos::HostSpace>& x) const
{
  double y;
  double xv = x[xi_];
  int n = x_.extent(0);
  if (xv <= x_.view_host()[0]) {
    y = y_.view_host()[0];
  } else if (xv > x_.view_host()[n - 1]) {
    y = y_.view_host()[n - 1];
  } else {
    // binary search to find interval containing xv
    int j1 = 0, j2 = n - 1;
    while (j2 - j1 > 1) {
      int j = (j1 + j2) / 2;
      // if (xv >= x_[j]) { // right continuous
      if (xv > x_.view_host()[j]) { // left continuous
        j1 = j;
      } else {
        j2 = j;
      }
    }
    // Now have x_[j1] <= xv < x_[j2], if right continuous
    // or x_[j1] < xv <= x_[j2], if left continuous
    switch (form_.view_host()[j1]) {
    case LINEAR:
      // Linear interpolation between x[j1] and x[j2]
      y = y_.view_host()[j1] + ((y_.view_host()[j2] - y_.view_host()[j1]) / (x_.view_host()[j2] - x_.view_host()[j1])) * (xv - x_.view_host()[j1]);
      break;
    case CONSTANT:
      y = y_.view_host()[j1];
      break;
    case FUNCTION:
      y = (*func_[j1])(x);
    }
  }
  return y;
}

} // namespace Amanzi
