/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Dylan Harp (dharp@lanl.gov)
*/

//! FunctionBilinear: a piecewise bilinear function.

/*!

A piecewise bilinear function extends the linear form of the tabular function to
two variables.

Define :math:`i(x) = i : x_i < x <= x_{{i+1}}` and similarly :math:`j(y) = j :
y_j < y <= y_{{j+1}}` for monotonically increasing :math:`x_i` and :math:`y_j`.

Given a two-dimensional array :math:`u_{{i,j}}`, :math:`f` is then defined by
bilinear interpolation on :math:`u_{{i(x),j(y)}}, u_{{i(x)+1,j(y)}},
u_{{i(x),j(y)+1}}, u_{{i(x)+1,j(y)+1}}, if :math:`(x,y)` is in
:math:`[x_0,x_n] \times [y_0,y_m]`, linear interpolation if one of :math:`x,y`
are out of those bounds, and constant at the corner value if both are out of
bounds.

* `"file`" ``[string]`` HDF5 filename of the data
* `"row header`" ``[string]`` name of the row dataset, the :math:`x_i`
* `"row coordinate`" ``[string]`` one of `"t`",`"x`",`"y`",`"z`"
* `"column header`" ``[string]`` name of the column dataset, the :math:`y_i`
* `"column coordinate`" ``[string]`` one of `"t`",`"x`",`"y`",`"z`"
* `"value header`" ``[string]`` name of the values dataset, the
:math:`u_{{i,j}}`

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
#ifndef AMANZI_BILINEAR_FUNCTION_HH_
#define AMANZI_BILINEAR_FUNCTION_HH_

#include <vector>

#include "Function.hh"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Amanzi {

class FunctionBilinear : public Function {
 public:
  FunctionBilinear(const Kokkos::View<double*,Kokkos::HostSpace>& x,
                   const Kokkos::View<double*,Kokkos::HostSpace>& y,
                   const Kokkos::View<double**,Kokkos::HostSpace>& v, const int xi, const int yi);
  ~FunctionBilinear(){};
  FunctionBilinear* Clone() const { return new FunctionBilinear(*this); }
  double operator()(const Kokkos::View<double*,Kokkos::HostSpace>&) const;

  KOKKOS_INLINE_FUNCTION double
  apply_gpu(const Kokkos::View<double**>& x, const int& i) const
  {
    double v;
    int nx = x_.extent(0);
    int ny = y_.extent(0);
    double xv = x(xi_, i);
    double yv = x(yi_, i);
    auto vv = v_.view_device(); 
    auto vx = x_.view_device(); 
    auto vy = y_.view_device(); 
    // if xv and yv are out of bounds
    if (xv <= vx[0] && yv <= vy[0]) {
      v = vv(0, 0);
    } else if (xv >= vx[nx - 1] && yv <= vy[0]) {
      v = vv(nx - 1, 0);
    } else if (xv >= vx[nx - 1] && yv >= vy[ny - 1]) {
      v = vv(nx - 1, ny - 1);
    } else if (xv <= vx[0] && yv >= vy[ny - 1]) {
      v = vv(0, ny - 1);
    } else {
      // binary search to find interval containing xv
      int j1 = 0, j2 = nx - 1;
      while (j2 - j1 > 1) {
        int j = (j1 + j2) / 2;
        if (xv >= vx[j]) { // right continuous
                           // if (xv > vx[j]) { // left continuous
          j1 = j;
        } else {
          j2 = j;
        }
      }
      // binary search to find interval containing yv
      int k1 = 0, k2 = ny - 1;
      while (k2 - k1 > 1) {
        int k = (k1 + k2) / 2;
        if (yv >= vy[k]) { // right continuous
                           // if (yv > vy[k]) { // left continuous
          k1 = k;
        } else {
          k2 = k;
        }
      }
      // if only xv is out of bounds, linear interpolation
      if (xv <= vx[0] && yv > vy[0] && yv < vy[ny - 1]) {
        v = vv(0, k1) +
            ((vv(0, k2) - vv(0, k1)) / (vy[k2] - vy[k1])) * (yv - vy[k1]);
      } else if (xv > vx[nx - 1] && yv > vy[0] && yv < vy[ny - 1]) {
        v = vv(nx - 1, k1) +
            ((vv(nx - 1, k2) - vv(nx - 1, k1)) / (vy[k2] - vy[k1])) *
              (yv - vy[k1]);
        // if only yv is out of bounds, linear interpolation
      } else if (yv <= vy[0] && xv > vx[0] && xv < vx[nx - 1]) {
        v = vv(j1, 0) +
            ((vv(j2, 0) - vv(j1, 0)) / (vx[j2] - vx[j1])) * (xv - vx[j1]);
      } else if (yv > vy[ny - 1] && xv > vx[0] && xv < vx[nx - 1]) {
        v = vv(j1, ny - 1) +
            ((vv(j2, ny - 1) - vv(j1, ny - 1)) / (vx[j2] - vx[j1])) *
              (xv - vx[j1]);
      } else {
        // bilinear interpolation
        v = vv(j1, k1) * (vx[j2] - xv) * (vy[k2] - yv) +
            vv(j2, k1) * (xv - vx[j1]) * (vy[k2] - yv) +
            vv(j1, k2) * (vx[j2] - xv) * (yv - vy[k1]) +
            vv(j2, k2) * (xv - vx[j1]) * (yv - vy[k1]);
        v = v / ((vx[j2] - vx[j1]) * (vy[k2] - vy[k1]));
      }
    }
    return v;
  }

  void apply(const Kokkos::View<double**>& in, Kokkos::View<double*>& out) const
  {
    assert(in.extent(1) == out.extent(0));
    Kokkos::parallel_for(
      "FunctionBilinear::apply",
      in.extent(1), KOKKOS_LAMBDA(const int& i) { out(i) = apply_gpu(in, i); });
  }

 private:
  Kokkos::DualView<double*> x_, y_;
  Kokkos::DualView<double**> v_;
  int xi_, yi_;

 private: // helper functions
  void check_args(const Kokkos::View<double*,Kokkos::HostSpace>&, const Kokkos::View<double*,Kokkos::HostSpace>&,
                  const Kokkos::View<double**,Kokkos::HostSpace>&) const;
};

} // namespace Amanzi

#endif // AMANZI_BILINEAR_FUNCTION_HH_
