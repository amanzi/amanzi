/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*! Input/Output via NetCDF file for extruded-column (MPAS-style) meshes,
  or their 2D surface meshes.

  For a 3D mesh, scoped to extruded/columnar meshes whose columns are
  uniform-height and memory-contiguous -- mesh->columns must exist
  (buildColumns() already called, e.g. via MeshFactory's "build columns"
  parameter), every column must have the same cell count, and each
  column's cells must be consecutive global ids.  ColumnarMeshInfo
  verifies all three directly against mesh->columns and errors clearly
  if any fails; see its .cc for details.

  For a 2D mesh (mesh->getManifoldDimension() == 2, e.g. a surface mesh --
  note this is the topological/manifold dimension, not getSpaceDimension(),
  since a surface mesh extracted from a 3D domain is still embedded in
  3-space), there is no vertical structure to ravel at all -- no columns
  are required, and CELL fields are written flat/unraveled, the same as
  any other location.  See ColumnarMeshInfo::is_columnar.

  All timesteps of a run are written as time slices into one on-disk file,
  named `"<file name base>.nc"` -- createTimestep/finalizeTimestep bracket
  a single write at the next time index, they do not open or close the
  file.

  write(attrs, int/double) is time-indexed once vis has begun: any such
  scalar written before the first createTimestep() call (e.g. a constant
  like "time units") is instead a plain, time-independent attribute.  See
  time_index_'s comment below.

  Per field, the DoF axis is either split into separate named/numbered
  variables (when attrs carries "subfieldnames") or kept as a genuine 4th
  dimension on one variable.

  Only CELL-location fields get the columnar reshape below; any other
  location (node, face, boundary_face, ...) is written flat/unraveled
  instead -- see write_'s comment.  The location is required: attrs must
  carry an AmanziMesh::Entity_kind "location" parameter (set by
  Output::writeCompositeVector_ from the CompositeVector's own
  per-component Entity_kind) -- write(Vector_type&)/write(MultiVector_type&)
  throw clearly if it is missing rather than guessing one.

  The ravel/unravel between flat mesh cell ordering and (column, layer) is
  closed-form arithmetic given the uniform-height/contiguity properties
  ColumnarMeshInfo verifies:
    column = global_cell_id / ncells_per_column
    layer  = global_cell_id % ncells_per_column
  and the column-axis I/O offset for hyperslabs is column_offset, an
  exclusive prefix sum of ncolumns_local across ranks.

.. _output-vertically-structured-netcdf-spec:
.. admonition:: output-vertically-structured-netcdf-spec

   * `"file name base`" ``[string]`` Base path for the on-disk file --
     written to `"<file name base>.nc"`.  All timesteps of the run are
     written to this one file.  OutputFactory::createForVis already
     appends the domain name to this parameter before construction, so
     this class itself does no domain-name handling.
*/

#pragma once

#include <string>
#include <memory>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "AmanziTypes.hh"
#include "Mesh.hh"
#include "Output.hh"
#include "Input.hh"
#include "FileNetCDF.hh"

namespace Amanzi {

namespace Impl {

// Shared columnar-mesh bookkeeping.  Two cases, keyed on
// mesh->getManifoldDimension() (topological dimension -- a surface mesh
// extracted from a 3D domain is still embedded in 3-space, so
// getSpaceDimension() would wrongly read 3 for it; manifold dimension is
// the right query):
//
//   3D (is_columnar true): mesh->columns must exist (buildColumns()
//   already called) and be uniform-height/memory-contiguous (the
//   properties FileNetCDF's columnar ravel depends on) -- verified
//   directly against mesh->columns, erroring clearly if either fails.
//   ncells_per_column/ncolumns_global/ncolumns_local/column_offset are
//   populated for the ravel/unravel math.
//
//   2D (is_columnar false, a surface mesh): no columns, no ravel --
//   CELL fields are written flat/unraveled instead, same as any other
//   location.  ncells_per_column/ncolumns_*/column_offset are left at
//   -1; callers must check is_columnar before using them.
//
// Any other manifold dimension is an error.  Shared by
// InputVerticallyStructuredNetCDF and OutputVerticallyStructuredNetCDF
// so the validation/arithmetic is written exactly once.
struct ColumnarMeshInfo {
  explicit ColumnarMeshInfo(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  Teuchos::RCP<const AmanziMesh::Mesh> mesh;
  bool is_columnar;
  LO ncells_per_column;
  GO ncolumns_global;
  GO column_offset; // this rank's first owned column (global index)
  LO ncolumns_local;
};

} // namespace Impl


// Only readRange() is implemented (used by FileTimeInterpolator to read
// standalone time-series input files, e.g. forcing data) -- every other
// Input method falls back to Input's default, which throws, since
// nothing currently reads a value through them for this file type.
class InputVerticallyStructuredNetCDF : public Input {
 public:
  InputVerticallyStructuredNetCDF(const Comm_ptr_type& comm,
                                  const std::string& filename,
                                  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  virtual void
  readRange(const Teuchos::ParameterList& attrs,
            std::size_t time_start,
            std::size_t time_count,
            Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> buf) const override;

 protected:
  std::unique_ptr<Utils::FileNetCDF> file_;
  Impl::ColumnarMeshInfo info_;
};


class OutputVerticallyStructuredNetCDF : public Output {
 public:
  OutputVerticallyStructuredNetCDF(Teuchos::ParameterList& plist,
                                   const Comm_ptr_type& comm,
                                   const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  virtual void createTimestep(double time, int cycle) override;
  virtual void finalizeTimestep() override;
  virtual std::string getFilename(int cycle) const override;

  virtual void write(const Teuchos::ParameterList& attrs, const int& val) const override;
  virtual void write(const Teuchos::ParameterList& attrs, const double& val) const override;
  virtual void write(const Teuchos::ParameterList& attrs, const std::string& val) const override;
  // both require attrs' "location" parameter -- see the class comment.
  virtual void write(const Teuchos::ParameterList& attrs, const Vector_type& vec) const override;
  virtual void write(const Teuchos::ParameterList& attrs, const IntVector_type& vec) const override;

  virtual void
  write(const Teuchos::ParameterList& attrs, const MultiVector_type& vec) const override;

 protected:
  // Writes vec at the current time index: columnar-reshaped ([ncolumn,
  // ncell_per_column], see the class comment) if location is CELL and
  // info_.is_columnar (a 3D, columnar mesh); otherwise written
  // flat/unraveled ([N], N = vec's own local entity count), the flat
  // dimension named "n" + location (e.g. "nboundary_face", "nface",
  // "ncell" for a CELL field on a 2D surface mesh -- matching
  // ncolumn/ncell_per_column's own "n"+singular-noun style, and UGRID's
  // convention for mesh-entity-count dimensions) so distinct fields at
  // the same location share one dimension -- node, face, boundary_face,
  // and CELL-on-a-surface-mesh have no per-column list in MeshColumns
  // and so cannot be reshaped into columns at all.
  void write_(const std::string& varname,
              const MultiVector_type& vec,
              AmanziMesh::Entity_kind location) const;
  void writeVerticallyStructured_(const std::string& varname, const MultiVector_type& vec) const;
  void writeFlat_(const std::string& varname,
                  const MultiVector_type& vec,
                  AmanziMesh::Entity_kind location) const;

  std::string filename_;
  Comm_ptr_type comm_;

  int cycle_;
  double time_;

  // -1 until the first createTimestep() call; a scalar written while
  // time_index_ < 0 is a time-independent attribute (e.g. "time units"),
  // written with FileNetCDF::writeAttribute(name, value).  Once vis has
  // begun, every scalar is a per-timestep value, written with
  // FileNetCDF::writeAttribute(name, value, time_index_).
  int time_index_;

  std::unique_ptr<Utils::FileNetCDF> file_;
  Impl::ColumnarMeshInfo info_;
};


} // namespace Amanzi
