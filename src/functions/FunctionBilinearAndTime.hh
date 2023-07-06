/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! FunctionNDLinear: a piecewise ND-linear function.
/*!

A piecewise bilinear function that is additionally linear interpolated in time.

Define :math:`i(x) = i : x_i < x <= x_{{i+1}}` and similarly :math:`j(y) = j : y_j < y <= y_{{j+1}}` 
for monotonically increasing :math:`x_i` and :math:`y_j`.

Given a two-dimensional array :math:`u_{{i,j}}`, :math:`f` is then defined by
bilinear interpolation on 
:math:`u_{i(x),j(y)}, u_{i(x)+1,j(y)}, u_{i(x),j(y)+1}, u_{i(x)+1,j(y)+1}`, if :math:`(x,y)` is in
:math:`[x_0,x_n] \times [y_0,y_m]`, linear interpolation if one of :math:`x,y`
are out of those bounds, and constant at the corner value if both are out of
bounds.

* `"file`" ``[string]`` HDF5 filename of the data
* `"time header`" ``[string]`` **time** Name of the temporal dimension indices, the :math:`t_i`.
* `"row header`" ``[string]`` **x** name of the row dataset, the :math:`x_i`
* `"row coordinate`" ``[string]`` **x** one of `"x`",`"y`",`"z`"
* `"column header`" ``[string]`` **y** name of the column dataset, the :math:`y_i`
* `"column coordinate`" ``[string]`` **y** one of `"x`",`"y`",`"z`"
* `"value header`" ``[string]`` name of the values dataset, the :math:`u_{{i,j}}`
* `"forms`" ``[string]`` **linear** Describes the temporal interpolant, one
  of `"linear`" or `"constant`", where `"linear`" is therefore trilinear
  interpolation (2x space and time) and `"constant`" indicates that the value
  on an interval is provided by the left point's (earlier in time) value.

Example1:

.. code-block:: xml

  <ParameterList name="function-nd-linear">
    <Parameter name="file" type="string" value="head.h5"/>
    <Parameter name="time header" type="string" value="time"/>
    <Parameter name="row header" type="string" value="x"/>
    <Parameter name="column header" type="string" value="y"/>
    <Parameter name="value header" type="string" value="rain"/>
  </ParameterList>

An example HDF5 file, called head.h5, might then look like:

|
| time: array(4) = [0, 60, 120, 180]
| x: array(3) = [0, 1, 2]
| y: array(2) = [0, 1]
| rain: (group)
|      | 0: array(3,2) = ...   # values at time 0
|      | 1: array(3,2) = ...   # values at time 60
|      | 2: array(3,2) = ...   # values at time 120
|      | 3: array(3,2) = ...   # values at time 180


*/
#ifndef AMANZI_BILINEAR_AND_TIME_FUNCTION_HH_
#define AMANZI_BILINEAR_AND_TIME_FUNCTION_HH_

#include <vector>

#include "FunctionBilinear.hh"

namespace Amanzi {

class FunctionBilinearAndTime : public Function {
 public:
  FunctionBilinearAndTime(const std::string& filename,
                          const std::string& time_header,
                          const std::string& row_header,
                          const std::string& row_coordinate,
                          const std::string& column_header,
                          const std::string& column_coordinate,
                          const std::string& val_header,
                          Form_kind form);

  FunctionBilinearAndTime(const FunctionBilinearAndTime& other);

  ~FunctionBilinearAndTime(){};
  std::unique_ptr<Function> Clone() const
  {
    return std::make_unique<FunctionBilinearAndTime>(*this);
  }

  double operator()(const std::vector<double>& x) const;

 private:
  std::unique_ptr<FunctionBilinear> Load_(const int i) const;

 private:
  std::string row_header_, col_header_, val_header_;
  int row_index_, col_index_;
  std::vector<double> times_;
  std::string filename_;
  Form_kind form_;

  mutable double t_before_;
  mutable double t_after_;
  mutable int current_interval_;

  mutable std::unique_ptr<Function> val_before_;
  mutable std::unique_ptr<Function> val_after_;
};

} // namespace Amanzi

#endif // AMANZI_BILINEAR_FUNCTION_HH_
