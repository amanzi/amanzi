/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Dylan Harp (dharp@lanl.gov)
*/

//! FunctionBilinear: a piecewise bilinear function.
/*!

A piecewise bilinear function extends the linear form of the tabular function to
two variables.

Define :math:`i(x) = i : x_i < x <= x_{{i+1}}` and similarly :math:`j(y) = j :
y_j < y <= y_{{j+1}}` for monotonically increasing :math:`x_i` and :math:`y_j`.

Given a two-dimensional array :math:`u_{i,j}`, :math:`f` is then defined by
bilinear interpolation on :math:`u_{i(x),j(y)}, u_{i(x)+1,j(y)}, u_{i(x),j(y)+1}, u_{i(x)+1,j(y)+1}`,
if :math:`(x,y)` is in :math:`[x_0,x_n] \times [y_0,y_m]`, linear interpolation if one of :math:`x,y`
are out of those bounds, and constant at the corner value if both are out of
bounds.

.. _function-bilinear-spec:
.. admonition:: function-bilinear-spec

   * `"file`" ``[string]`` HDF5 filename of the data
   * `"row header`" ``[string]`` name of the row dataset, the :math:`x_i`
   * `"row coordinate`" ``[string]`` one of `"t`",`"x`",`"y`",`"z`"
   * `"column header`" ``[string]`` name of the column dataset, the :math:`y_i`
   * `"column coordinate`" ``[string]`` one of `"t`",`"x`",`"y`",`"z`"
   * `"value header`" ``[string]`` name of the values dataset, the :math:`u_{{i,j}}`

Example:

.. code-block:: xml

  <ParameterList name="function-bilinear">
    <Parameter name="file" type="string" value="pressure.h5"/>
    <Parameter name="row header" type="string" value="/time"/>
    <Parameter name="row coordinate" type="string" value="t"/>
    <Parameter name="column header" type="string" value="/x"/>
    <Parameter name="column coordinate" type="string" value="x"/>
    <Parameter name="value header" type="string" value="/pressure"/>
  </ParameterList>

*/
#pragma once

#include <vector>

#include "Function.hh"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Amanzi {


class FunctionBilinear : public Function {
 public:
  FunctionBilinear(const Kokkos::View<const double*, Kokkos::HostSpace>& x,
                   const Kokkos::View<const double*, Kokkos::HostSpace>& y,
                   const Kokkos::View<const double**, Kokkos::HostSpace>& v,
                   const int xi,
                   const int yi);

  std::unique_ptr<Function> Clone() const override { return std::make_unique<FunctionBilinear>(*this); }

  double operator()(const Kokkos::View<const double**, Kokkos::HostSpace>& in) const override;

  void apply(const Kokkos::View<const double**>& in,
             Kokkos::View<double*>& out,
             const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const override;

 private:
  Kokkos::DualView<const double*> x_, y_;
  Kokkos::DualView<const double**> v_;
  int xi_, yi_;

 private: // helper functions
  void checkArgs_(const Kokkos::View<const double*, Kokkos::HostSpace>&,
                  const Kokkos::View<const double*, Kokkos::HostSpace>&,
                  const Kokkos::View<const double**, Kokkos::HostSpace>&) const;
};


namespace Impl {

template<class XYView_type,
         class VView_type,
         class InView_type>
class FunctionBilinearFunctor {
 public:
  FunctionBilinearFunctor(const XYView_type& x,
                          const XYView_type& y,
                          size_t xi,
                          size_t yi,
                          const VView_type& v,
                          const InView_type& in)
    : x_(x), y_(y), xi_(xi), yi_(yi), v_(v), in_(in) {}

  KOKKOS_INLINE_FUNCTION
  double operator()(const int i) const {
    int nx = x_.extent(0);
    int ny = y_.extent(0);
    double xv = in_(xi_, i);
    double yv = in_(yi_, i);

    // if xv and yv are out of bounds
    double v = 0.;
    if (xv <= x_[0] && yv <= y_[0]) {
      v = v_(0, 0);
    } else if (xv >= x_[nx - 1] && yv <= y_[0]) {
      v = v_(nx - 1, 0);
    } else if (xv >= x_[nx - 1] && yv >= y_[ny - 1]) {
      v = v_(nx - 1, ny - 1);
    } else if (xv <= x_[0] && yv >= y_[ny - 1]) {
      v = v_(0, ny - 1);
    } else {
      // find interval containing xv, yv
      int j2 = 0;
      while ((j2 < nx) && (xv > x_[j2])) ++j2;
      int j1 = j2 - 1;

      int k2 = 0;
      while ((k2 < ny) && (yv > y_[k2])) ++k2;
      int k1 = k2 - 1;

      // if only xv is out of bounds, linear interpolation
      if (j2 == 0) {
        v = v_(j2, k1) + ((v_(j2, k2) - v_(j2, k1)) / (y_[k2] - y_[k1])) * (yv - y_[k1]);
      } else if (j2 == nx) {
        v = v_(j1, k1) + ((v_(j1, k2) - v_(j1, k1)) / (y_[k2] - y_[k1])) * (yv - y_[k1]);
        // if only yv is out of bounds, linear interpolation
      } else if (k2 == 0) {
        v = v_(j1, k2) + ((v_(j2, k2) - v_(j1, k2)) / (x_[j2] - x_[j1])) * (xv - x_[j1]);
      } else if (k2 == ny) {
        v = v_(j1, k1) + ((v_(j2, k1) - v_(j1, k1)) / (x_[j2] - x_[j1])) * (xv - x_[j1]);
      } else {
        // bilinear interpolation
        v = v_(j1, k1) * (x_[j2] - xv) * (y_[k2] - yv) +
            v_(j2, k1) * (xv - x_[j1]) * (y_[k2] - yv) +
            v_(j1, k2) * (x_[j2] - xv) * (yv - y_[k1]) + v_(j2, k2) * (xv - x_[j1]) * (yv - y_[k1]);
        v = v / ((x_[j2] - x_[j1]) * (y_[k2] - y_[k1]));
      }
    }
    return v;
  }

 private:
  XYView_type x_, y_;
  int xi_, yi_;
  VView_type v_;
  InView_type in_;
};

} // namespace Impl
} // namespace Amanzi

