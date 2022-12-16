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

The functional forms :math:`{f_j}` may be constant, which uses the left endpoint, i.e.

:math:`f_i(x) = y_i`,

linear, i.e.

:math:`f_i(x) = ( y_i * (x - x_i) + y_{{i+1}} * (x_{{i+1}} - x) ) / (x_{{i+1}} - x_i)`

or arbitrary, in which the :math:`f_j` must be provided.

The :math:`x_i` and :math:`y_i` may be provided in one of two ways -- explicitly in the input spec or from an HDF5 file.  The length of these must be equal, and the :math:`x_i` must be monotonically increasing.  Forms, as defined on intervals, must be of length equal to the length of the :math:`x_i` less one.

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
    <Parameter name="x values" type="Array(double)" value="{0.0, 1.0, 2.0, 3.0}"/>
    <Parameter name="x coordinate" type="string" value="t"/>
    <Parameter name="y values" type="Array(double)" value="{0.0, 1.0, 2.0, 2.0}"/>
    <Parameter name="forms" type="Array(string)" value="{linear, constant, USER_FUNC}"/>

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

The example below would perform linear-interpolation on the intervals provided by data within the hdf5 file `"my_data.h5`".

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

namespace Amanzi {

class FunctionTabular : public Function {
 public:
  enum Form { LINEAR, CONSTANT, FUNCTION };

 public:
  FunctionTabular(const std::vector<double>& x, const std::vector<double>& y, const int xi);
  FunctionTabular(const std::vector<double>& x,
                  const std::vector<double>& y,
                  const int xi,
                  const std::vector<Form>& form);
  FunctionTabular(const std::vector<double>& x,
                  const std::vector<double>& y,
                  const int xi,
                  const std::vector<Form>& form,
                  std::vector<std::unique_ptr<Function>> func);
  FunctionTabular(const FunctionTabular& other);
  ~FunctionTabular(){};

  std::unique_ptr<Function> Clone() const { return std::make_unique<FunctionTabular>(*this); }
  double operator()(const std::vector<double>& x) const;

 private:
  std::vector<double> x_, y_;
  int xi_;
  std::vector<Form> form_;
  std::vector<std::unique_ptr<Function>> func_;

 private: // helper functions
  void check_args(const std::vector<double>&,
                  const std::vector<double>&,
                  const std::vector<Form>&) const;
};

} // namespace Amanzi

#endif // AMANZI_TABULAR_FUNCTION_HH_
