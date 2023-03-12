/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Dylan Harp (dharp@lanl.gov)
*/

//! FunctionBilinear: a piecewise bilinear function.
/*!

A piecewise bilinear function extends the linear form of the tabular function to two variables.

Define :math:`i(x) = i : x_i < x <= x_{{i+1}}` and similarly :math:`j(y) = j : y_j < y <= y_{{j+1}}` for monotonically increasing :math:`x_i` and :math:`y_j`.

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
#ifndef AMANZI_BILINEAR_FUNCTION_HH_
#define AMANZI_BILINEAR_FUNCTION_HH_

#include <vector>

#include "Function.hh"
#include "Epetra_SerialDenseMatrix.h"

namespace Amanzi {

class FunctionBilinear : public Function {
 public:
  FunctionBilinear(const std::vector<double>& x,
                   const std::vector<double>& y,
                   const Epetra_SerialDenseMatrix& v,
                   const int xi,
                   const int yi);
  ~FunctionBilinear(){};
  std::unique_ptr<Function> Clone() const { return std::make_unique<FunctionBilinear>(*this); }
  double operator()(const std::vector<double>& x) const;

 private:
  std::vector<double> x_, y_;
  Epetra_SerialDenseMatrix v_;
  int xi_, yi_;

 private: // helper functions
  void check_args(const std::vector<double>&,
                  const std::vector<double>&,
                  const Epetra_SerialDenseMatrix&) const;
};

} // namespace Amanzi

#endif // AMANZI_BILINEAR_FUNCTION_HH_
