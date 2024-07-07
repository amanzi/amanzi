/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//!

#include "errors.hh"
#include "ViewUtils.hh"

#include "FunctionTabular.hh"

namespace Amanzi {

FunctionTabular::FunctionTabular(const Kokkos::View<const double*, Kokkos::HostSpace>& x,
                                 const Kokkos::View<const double*, Kokkos::HostSpace>& y,
                                 const int xi)
  : xi_(xi)
{
  Kokkos::View<Form_kind*, Kokkos::HostSpace> form("form", x.size()-1);
  Kokkos::deep_copy(form, Form_kind::LINEAR);
  check_args(x, y, form);

  x_ = asDualView(x);
  y_ = asDualView(y);
  form_ = asDualView(form);
}

FunctionTabular::FunctionTabular(const Kokkos::View<const double*, Kokkos::HostSpace>& x,
                                 const Kokkos::View<const double*, Kokkos::HostSpace>& y,
                                 const int xi,
                                 const Kokkos::View<const Form_kind*, Kokkos::HostSpace>& form)
  : xi_(xi)
{
  check_args(x, y, form);
  x_ = asDualView(x);
  y_ = asDualView(y);
  form_ = asDualView(form);
}


void
FunctionTabular::check_args(const Kokkos::View<const double*, Kokkos::HostSpace>& x,
                            const Kokkos::View<const double*, Kokkos::HostSpace>& y,
                            const Kokkos::View<const Form_kind*, Kokkos::HostSpace>& form) const
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


void
FunctionTabular::apply(const Kokkos::View<const double**>& in,
                       Kokkos::View<double*>& out,
                       const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const
{
  auto f = Impl::FunctionTabularFunctor(x_.view_device(), y_.view_device(), form_.view_device(), xi_, in);

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

double
FunctionTabular::operator()(const Kokkos::View<const double**, Kokkos::HostSpace>& x) const
{
  auto f = Impl::FunctionTabularFunctor(x_.view_device(), y_.view_device(), form_.view_device(), xi_, x);
  return f(0);
}

} // namespace Amanzi
