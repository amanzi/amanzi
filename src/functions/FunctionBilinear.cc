/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//!
#include "errors.hh"
#include "FunctionBilinear.hh"
#include "ViewUtils.hh"

namespace Amanzi {

FunctionBilinear::FunctionBilinear(const Kokkos::View<const double*, Kokkos::HostSpace>& x,
                                   const Kokkos::View<const double*, Kokkos::HostSpace>& y,
                                   const Kokkos::View<const double**, Kokkos::HostSpace>& v,
                                   const int xi,
                                   const int yi)
  : xi_(xi), yi_(yi)
{
  checkArgs_(x, y, v);
  x_ = asDualView(x);
  y_ = asDualView(y);
  v_ = asDualView(v);
}


void
FunctionBilinear::checkArgs_(const Kokkos::View<const double*, Kokkos::HostSpace>& x,
                             const Kokkos::View<const double*, Kokkos::HostSpace>& y,
                             const Kokkos::View<const double**, Kokkos::HostSpace>& v) const
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


double
FunctionBilinear::operator()(const Kokkos::View<const double**, Kokkos::HostSpace>& in) const {
  auto f = Impl::FunctionBilinearFunctor(x_.view_host(), y_.view_host(), xi_, yi_, v_.view_host(), in);
  return f(0);
}


void
FunctionBilinear::apply(const Kokkos::View<const double**>& in,
           Kokkos::View<double*>& out,
           const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const
{
  auto f = Impl::FunctionBilinearFunctor(x_.view_device(), y_.view_device(), xi_, yi_, v_.view_device(), in);

  if (ids) {
    auto ids_loc = *ids;
    Kokkos::parallel_for(
      "FunctionBilinear::apply1", ids_loc.extent(0), KOKKOS_LAMBDA(const int i) {
        out(ids_loc(i)) = f(ids_loc(i));
      });
  } else {
    Kokkos::parallel_for("FunctionBilinear::apply2", in.extent(1), KOKKOS_LAMBDA(const int i) {
        out(i) = f(i);
      });
    }
  }

} // namespace Amanzi
