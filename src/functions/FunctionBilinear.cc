/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//!
#include "FunctionBilinear.hh"
#include "errors.hh"

namespace Amanzi {

FunctionBilinear::FunctionBilinear(const Kokkos::View<double*, Kokkos::HostSpace>& x,
                                   const Kokkos::View<double*, Kokkos::HostSpace>& y,
                                   const Kokkos::View<double**, Kokkos::HostSpace>& v,
                                   const int xi,
                                   const int yi)
  : xi_(xi), yi_(yi)
{
  checkArgs_(x, y, v);

  Kokkos::resize(x_, x.extent(0));
  Kokkos::resize(y_, y.extent(0));
  Kokkos::resize(v_, v.extent(0), v.extent(1));

  Kokkos::deep_copy(x_.view_host(), x);
  Kokkos::deep_copy(y_.view_host(), y);
  Kokkos::deep_copy(v_.view_host(), v);

  Kokkos::deep_copy(x_.view_device(), x_.view_host());
  Kokkos::deep_copy(y_.view_device(), y_.view_host());
  Kokkos::deep_copy(v_.view_device(), v_.view_host());

}


void
FunctionBilinear::checkArgs_(const Kokkos::View<double*, Kokkos::HostSpace>& x,
                             const Kokkos::View<double*, Kokkos::HostSpace>& y,
                             const Kokkos::View<double**, Kokkos::HostSpace>& v) const
{
  if (x.extent(0) != v.extent(0)) {
    Errors::Message m;
    m << "the number of x values and row in v differ";
    Exceptions::amanzi_throw(m);
  }
  if (y.extent(0) != v.extent(1)) {
    Errors::Message m;
    m << "the number of y values and columns in v differ";
    Exceptions::amanzi_throw(m);
  }
  if (x.extent(0) < 2) {
    Errors::Message m;
    m << "at least two rows values must be given";
    Exceptions::amanzi_throw(m);
  }
  if (y.extent(0) < 2) {
    Errors::Message m;
    m << "at least two column values must be given";
    Exceptions::amanzi_throw(m);
  }
  for (int j = 1; j < x.extent(0); ++j) {
    if (x[j] <= x[j - 1]) {
      Errors::Message m;
      m << "x values are not strictly increasing";
      Exceptions::amanzi_throw(m);
    }
  }
  for (int j = 1; j < y.extent(0); ++j) {
    if (y[j] <= y[j - 1]) {
      Errors::Message m;
      m << "y values are not strictly increasing";
      Exceptions::amanzi_throw(m);
    }
  }
}

} // namespace Amanzi
