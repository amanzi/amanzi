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

.. _function-bilinear-and-time-spec:
.. admonition:: function-bilinear-and-time-spec

   * `"file`" ``[string]`` filename of the data
   * `"file format`" ``[string]`` **optional** `"HDF5`" or `"NetCDF`" --
     inferred from the filename's extension (".h5" or ".nc") if not provided
   * `"time header`" ``[string]`` **time** Name of the temporal dimension indices, the :math:`t_i`.
   * `"row header`" ``[string]`` **x** name of the row dataset, the :math:`x_i`
   * `"row coordinate`" ``[string]`` **x** one of `"x`",`"y`",`"z`"
   * `"column header`" ``[string]`` **y** name of the column dataset, the :math:`y_i`
   * `"column coordinate`" ``[string]`` **y** one of `"x`",`"y`",`"z`"
   * `"value header`" ``[string]`` name of the values dataset, the :math:`u_{{i,j}}`
   * `"forms`" ``[string]`` **linear** Describes the temporal interpolant, one
     of `"linear`" or `"constant`", where `"linear`" is therefore trilinear
     interpolation (2x space and time) and `"constant`" indicates that the
     value on an interval is provided by the left point's (earlier in time)
     value.
   * `"page size`" ``[int]`` **optional, default=2** Number of consecutive
     time slices kept resident at once.

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
| rain: array(4,3,2) -- one [3,2] slice per time index, in that order.

*/
#ifndef AMANZI_BILINEAR_AND_TIME_FUNCTION_HH_
#define AMANZI_BILINEAR_AND_TIME_FUNCTION_HH_

#include <memory>
#include <string>
#include <vector>

#include "AmanziComm.hh"
#include "AmanziTypes.hh"
#include "IOFileFormat.hh"
#include "FileReader.hh"
#include "TimeSeriesInterpolator.hh"

#include "FunctionBilinear.hh"

namespace Amanzi {

class FunctionBilinearAndTime : public Function {
 public:
  using Page_type = Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace>;

  FunctionBilinearAndTime(const std::string& filename,
                          const std::string& time_header,
                          const std::string& row_header,
                          const std::string& row_coordinate,
                          const std::string& column_header,
                          const std::string& column_coordinate,
                          const std::string& val_header,
                          Form_kind form,
                          int page_size = 2,
                          Utils::IOFileFormat_kind file_format = Utils::IOFileFormat_kind::HDF5,
                          const Comm_ptr_type& comm = getDefaultComm());

  FunctionBilinearAndTime(const FunctionBilinearAndTime& other);

  std::unique_ptr<Function> Clone() const override
  { return std::make_unique<FunctionBilinearAndTime>(*this); }

  double operator()(const Kokkos::View<const double**, Kokkos::HostSpace>& x) const override;

  void apply(
    const Kokkos::View<const double**>& in,
    Kokkos::View<double*>& out,
    const Kokkos::MeshView<const int*, Amanzi::DefaultMemorySpace>* ids = nullptr) const override;

 private:
  // (Re-)builds file_, row_/col_, and interp_ from the stored construction
  // parameters -- shared by the constructor and the copy constructor
  // (which cannot cheaply copy an open FileReader or a
  // TimeSeriesInterpolator whose load_page closure captures this).
  void initialize_();

  std::string filename_, time_header_, row_header_, col_header_, val_header_;
  std::string row_coordinate_, col_coordinate_;
  int row_index_, col_index_;
  Form_kind form_;
  int page_size_;
  Utils::IOFileFormat_kind file_format_;
  Comm_ptr_type comm_;

  std::unique_ptr<Utils::FileReader> file_;
  Kokkos::DualView<double*> row_, col_;
  int nx_, ny_;

  mutable std::unique_ptr<TimeSeriesInterpolator<Page_type>> interp_;
  mutable Page_type page_, result_;
  mutable Kokkos::View<double**, Kokkos::LayoutRight> result_d_;
};

} // namespace Amanzi

#endif // AMANZI_BILINEAR_FUNCTION_HH_
