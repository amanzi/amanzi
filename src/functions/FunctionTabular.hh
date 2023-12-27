/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//! FunctionTabular: Piecewise-defined function.
/*!

A piecewise function of one variable.

A tabular function is tabulated on a series of intervals; given values
:math:`{{x_i}}, {{y_i}},, i=0, ... n-1` and functional forms :math:`{{f_j}},,
j=0, ... n-2` a tabular function :math:`f(x)` is defined as:

.. math::
  \begin{matrix}
  f(x) &=& y_0, & x \le x_0,\\
  f(x) &=& f_{{i-1}}(x)  & x \in (x_{{i-1}}, x_i],\\
  f(x) &=& y_{{n-1}}, & x > x_{{n-1}}.
  \end{matrix}

The functional forms :math:`{f_j}` may be constant, which uses the left
endpoint, i.e.

:math:`f_i(x) = y_i`,

linear, i.e.

:math:`f_i(x) = ( y_i * (x - x_i) + y_{{i+1}} * (x_{{i+1}} - x) ) / (x_{{i+1}} -
x_i)`

or arbitrary, in which the :math:`f_j` must be provided.

The :math:`x_i` and :math:`y_i` may be provided in one of two ways -- explicitly
in the input spec or from an HDF5 file.  The length of these must be equal, and
the :math:`x_i` must be monotonically increasing.  Forms, as defined on
intervals, must be of length equal to the length of the :math:`x_i` less one.

Explicitly specifying the data:

.. _function-tabular-spec:
.. admonition:: function-tabular-spec

   * `"x values`" ``[Array(double)]`` the :math:`x_i`
   * `"y values`" ``[Array(double)]`` the :math:`y_i`

   * `"forms`" ``[Array(string)]`` **optional** Form of the interpolant, either
     `"constant`", `"linear`", or `"USER_DEFINED`" Default is linear for each *
     interval.  Note the length of this array should be one per interval, or
     one less than len of x and y values.

   * `"USER_DEFINED`" ``[function-spec]`` **optional** user-provided functional
     forms on the interval

   * `"x coordinate`" ``[string]`` **t**, `"x`", `"y`", `"z`" defines which
     coordinate direction the :math:`x_i` are formed, defaulting to time.

The below example defines a function that is zero on interval :math:`(-\infty,\,0]`,
linear on interval :math:`(0,\,1]`, constant (`f(x)=1`) on interval :math:`(1,\,2]`,
square root of `t` on interval :math:`(2,\,3]`,
and constant (`f(x)=2`) on interval :math:`(3,\,\infty]`.

Example:

.. code-block:: xml

  <ParameterList name="function-tabular">
    <Parameter name="x values" type="Array(double)"
value="{0.0, 1.0, 2.0, 3.0}"/> <Parameter name="x coordinate" type="string"
value="t"/> <Parameter name="y values" type="Array(double)"
value="{0.0, 1.0, 2.0, 2.0}"/> <Parameter name="forms" type="Array(string)"
value="{linear, constant, USER_FUNC}"/>

    <ParameterList name="USER_FUNC">
      <ParameterList name="function-standard-math">
        <Parameter name="operator" type="string" value="sqrt"/>
      </ParameterList>
    </ParameterList>
  </ParameterList>


Loading table from file.  (Note that `"USER_DEFINED`" is not an option here,
but could be made so if requested).

.. _function-tabular-fromfile-spec:
.. admonition:: function-tabular-fromfile-spec

   * `"file`" ``[string]`` filename of the HDF5 data
   * `"x header`" ``[string]`` name of the dataset for the :math:`x_i` in the file
   * `"y header`" ``[string]`` name of the dataset for the :math:`y_i` in the file
   * `"forms`" ``[string]`` **optional**, Form of the interpolant, either
     `"constant`" or `"linear`"

The example below would perform linear-interpolation on the intervals provided
by data within the hdf5 file `"my_data.h5`".

Example:

.. code-block:: xml

  <ParameterList name="function-tabular">
    <Parameter name="file" type="string" value="my_data.h5"/>
    <Parameter name="x coordinate" type="string" value="t"/>
    <Parameter name="x header" type="string" value="/time"/>
    <Parameter name="y header" type="string" value="/data"/>
  </ParameterList>

*/

#ifndef AMANZI_TABULAR_FUNCTION_HH_
#define AMANZI_TABULAR_FUNCTION_HH_

#include <memory>
#include <vector>

#include "Function.hh"
#include "AmanziTypes.hh"

namespace Amanzi {

class FunctionTabular : public Function {
 public:
  FunctionTabular(const Kokkos::View<double*, Kokkos::HostSpace>& x,
                  const Kokkos::View<double*, Kokkos::HostSpace>& y,
                  const int xi);
  FunctionTabular(const Kokkos::View<double*, Kokkos::HostSpace>& x,
                  const Kokkos::View<double*, Kokkos::HostSpace>& y,
                  const int xi,
                  const Kokkos::View<Form_kind*, Kokkos::HostSpace>& form);
  FunctionTabular(const Kokkos::View<double*, Kokkos::HostSpace>& x,
                  const Kokkos::View<double*, Kokkos::HostSpace>& y,
                  const int xi,
                  const Kokkos::View<Form_kind*, Kokkos::HostSpace>& form,
                  std::vector<std::unique_ptr<Function>> func);
  FunctionTabular(const FunctionTabular& other)
    : x_(other.x_), y_(other.y_), form_(other.form_), xi_(other.xi_)
  {
    for (const auto& f : other.func_) func_.emplace_back(f->Clone());
  }

  ~FunctionTabular(){};
  std::unique_ptr<Function> Clone() const { return std::make_unique<FunctionTabular>(*this); }
  double operator()(const Kokkos::View<double*, Kokkos::HostSpace>&) const;

  KOKKOS_INLINE_FUNCTION double apply_gpu(const Kokkos::View<double**>& x, const int i) const
  {
    double y;
    double xv = x(xi_, i);
    int n = x_.extent(0);
    if (xv <= x_.view_device()[0]) {
      y = y_.view_device()[0];
    } else if (xv > x_.view_device()[n - 1]) {
      y = y_.view_device()[n - 1];
    } else {
      // binary search to find interval containing xv
      int j1 = 0, j2 = n - 1;
      while (j2 - j1 > 1) {
        int j = (j1 + j2) / 2;
        // if (xv >= x_[j]) { // right continuous
        if (xv > x_.view_device()[j]) { // left continuous
          j1 = j;
        } else {
          j2 = j;
        }
      }
      // Now have x_[j1] <= xv < x_[j2], if right continuous
      // or x_[j1] < xv <= x_[j2], if left continuous
      switch (form_.view_device()[j1]) {
      case Form_kind::LINEAR:
        // Linear interpolation between x[j1] and x[j2]
        y = y_.view_device()[j1] + ((y_.view_device()[j2] - y_.view_device()[j1]) /
                                    (x_.view_device()[j2] - x_.view_device()[j1])) *
                                     (xv - x_.view_device()[j1]);
        break;
      case Form_kind::CONSTANT:
        y = y_.view_device()[j1];
        break;
      case Form_kind::FUNCTION:
        assert(false && "Not implemented for FUNCTION");
        //  y = (*func_[j1])(x);
      }
    }
    return y;
  }


  void apply(const Kokkos::View<double**>& in,
             Kokkos::View<double*>& out,
             const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const
  {
    if (ids) {
      auto ids_loc = *ids;
      Kokkos::parallel_for(
        "FunctionTabular::apply1", in.extent(1), KOKKOS_LAMBDA(const int& i) {
          out(ids_loc(i)) = apply_gpu(in, i);
        });
    } else {
      assert(in.extent(1) == out.extent(0));
      Kokkos::parallel_for(
        "FunctionTabular::apply2", in.extent(1), KOKKOS_LAMBDA(const int& i) {
          out(i) = apply_gpu(in, i);
        });
    }
  }

 private:
  Kokkos::DualView<double*, Amanzi::DeviceOnlyMemorySpace> x_;
  Kokkos::DualView<double*, Amanzi::DeviceOnlyMemorySpace> y_;
  Kokkos::DualView<Form_kind*, Amanzi::DeviceOnlyMemorySpace> form_;
  std::vector<std::unique_ptr<Function>> func_;
  int xi_;

 private: // helper functions
  void check_args(const Kokkos::View<double*, Kokkos::HostSpace>&,
                  const Kokkos::View<double*, Kokkos::HostSpace>&,
                  const Kokkos::View<Form_kind*, Kokkos::HostSpace>&) const;
};

} // namespace Amanzi

#endif // AMANZI_TABULAR_FUNCTION_HH_
