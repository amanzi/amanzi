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

#pragma once

#include <memory>
#include <vector>

#include "Function.hh"
#include "AmanziTypes.hh"

namespace Amanzi {

class FunctionTabular : public Function {
 public:
  FunctionTabular(const Kokkos::View<const double*, Kokkos::HostSpace>& x,
                  const Kokkos::View<const double*, Kokkos::HostSpace>& y,
                  const int xi);

  FunctionTabular(const Kokkos::View<const double*, Kokkos::HostSpace>& x,
                  const Kokkos::View<const double*, Kokkos::HostSpace>& y,
                  const int xi,
                  const Kokkos::View<const Form_kind*, Kokkos::HostSpace>& form);

  FunctionTabular(const FunctionTabular& other)
    : x_(other.x_), y_(other.y_), form_(other.form_), xi_(other.xi_)
  {}

  std::unique_ptr<Function> Clone() const override { return std::make_unique<FunctionTabular>(*this); }

  double operator()(const Kokkos::View<const double**, Kokkos::HostSpace>&) const override;

  void apply(const Kokkos::View<const double**>& in,
             Kokkos::View<double*>& out,
             const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids) const override;

 private:
  Kokkos::DualView<const double*> x_;
  Kokkos::DualView<const double*> y_;
  Kokkos::DualView<const Form_kind*> form_;
  int xi_;

 private: // helper functions
  void check_args(const Kokkos::View<const double*, Kokkos::HostSpace>&,
                  const Kokkos::View<const double*, Kokkos::HostSpace>&,
                  const Kokkos::View<const Form_kind*, Kokkos::HostSpace>&) const;
};


namespace Impl {

template <class DoubleView_type,
          class FormView_type,
          class InView_type>
class FunctionTabularFunctor {
 public:
  FunctionTabularFunctor(const DoubleView_type& x,
                        const DoubleView_type& y,
                        const FormView_type& form,
                        int xi,
                        const InView_type& in)
    : x_(x), y_(y), form_(form), xi_(xi), in_(in) {}

  KOKKOS_INLINE_FUNCTION
  double operator()(const int i) const
  {
    assert(xi_ >= 0 && xi_ < in_.extent(0));
    assert(i >= 0 && i < in_.extent(1));
    double xv = in_(xi_, i);
    int nx = x_.extent(0);

    double y(0.);
    if (xv <= x_[0]) {
      y = y_[0];
    } else if (xv > x_[nx - 1]) {
      y = y_[nx - 1];
    } else {
      int j2 = 0;
      while ((j2 < nx) && (xv > x_[j2])) ++j2;
      int j1 = j2 - 1;

      // Now have x_[j1] <= xv < x_[j2], if right continuous
      // or x_[j1] < xv <= x_[j2], if left continuous
      switch (form_[j1]) {
      case Form_kind::LINEAR:
        // Linear interpolation between x[j1] and x[j2]
        y = y_[j1] + ((y_[j2] - y_[j1]) /
                                    (x_[j2] - x_[j1])) *
                                     (xv - x_[j1]);
        break;
      case Form_kind::CONSTANT:
        y = y_[j1];
        break;
      case Form_kind::FUNCTION:
        assert(false && "Not implemented for FUNCTION");
      }
    }
    return y;
  }


 private:
  DoubleView_type x_, y_;
  FormView_type form_;
  int xi_;
  InView_type in_;
};

} // namespace Impl
} // namespace Amanzi


