/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Input/Output via NetCDF file for extruded-column (MPAS-style) meshes.

#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include "OutputUtils.hh"
#include "InputOutputVerticallyStructuredNetCDF.hh"

namespace Amanzi {

namespace InputOutputVerticallyStructuredNetCDF_Helpers {

ColumnarMeshInfo::ColumnarMeshInfo(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_)
  : mesh(mesh_), ncells_per_column(-1), ncolumns_global(-1), column_offset(-1), ncolumns_local(-1)
{
  // manifold (topological) dimension, not space (ambient) dimension -- a
  // surface mesh extracted from a 3D domain is still embedded in 3-space
  // (getSpaceDimension() == 3, its nodes carry x,y,z) but is topologically
  // 2D and has no vertical structure to ravel.
  int manifold_dim = mesh->getManifoldDimension();
  if (manifold_dim == 2) {
    // a surface mesh -- no vertical structure, no columns required; CELL
    // fields are written flat/unraveled, same as any other location.
    is_columnar = false;
    return;
  }
  if (manifold_dim != 3) {
    Errors::Message msg;
    msg << "InputOutputVerticallyStructuredNetCDF: mesh has manifold dimension " << manifold_dim
        << ", only 2 (surface) or 3 (columnar) are supported.";
    Exceptions::amanzi_throw(msg);
  }
  is_columnar = true;

  if (mesh->columns == Teuchos::null) {
    Errors::Message msg(
      "InputOutputVerticallyStructuredNetCDF: mesh->columns does not exist -- "
      "buildColumns() must be called (e.g. via MeshFactory's \"build columns\" parameter) "
      "before this mesh can be used for columnar NetCDF I/O.");
    Exceptions::amanzi_throw(msg);
  }

  ncolumns_local = mesh->columns->num_columns_owned;
  ncells_per_column =
    ncolumns_local > 0 ? (int)mesh->columns->getCells<MemSpace_kind::HOST>(0).extent(0) : 0;

  for (int col = 0; col != ncolumns_local; ++col) {
    auto col_cells = mesh->columns->getCells<MemSpace_kind::HOST>(col);
    if ((int)col_cells.extent(0) != ncells_per_column) {
      Errors::Message msg;
      msg << "InputOutputVerticallyStructuredNetCDF: column " << col << " has "
          << col_cells.extent(0) << " cells, column 0 has " << ncells_per_column
          << " -- columns must all be the same height for the columnar NetCDF ravel.";
      Exceptions::amanzi_throw(msg);
    }
    for (int i = 0; i != ncells_per_column - 1; ++i) {
      if (col_cells(i + 1) != col_cells(i) + 1) {
        Errors::Message msg;
        msg << "InputOutputVerticallyStructuredNetCDF: column " << col << "'s cells are not "
            << "memory-contiguous (cell " << i << " is local id " << col_cells(i) << ", cell "
            << (i + 1) << " is " << col_cells(i + 1)
            << ") -- required for the columnar NetCDF ravel.";
        Exceptions::amanzi_throw(msg);
      }
    }
  }

  long long local = ncolumns_local;
  long long global = 0, inclusive_scan = 0;
  Teuchos::reduceAll(*mesh->getComm(), Teuchos::REDUCE_SUM, 1, &local, &global);
  Teuchos::scan(*mesh->getComm(), Teuchos::REDUCE_SUM, 1, &local, &inclusive_scan);
  ncolumns_global = (int)global;
  column_offset = (int)(inclusive_scan - local);

  // Two independent, reduction-free cross-checks that ncells_per_column
  // (derived above from one rank's own column 0) is actually consistent
  // mesh-wide, not just self-consistent on this rank: this rank's owned
  // cells must divide evenly among its owned columns, and the mesh's
  // global cells must divide evenly among all columns, both landing on
  // the same ncells_per_column this rank already computed.
  if (ncolumns_local > 0) {
    int ncells_owned =
      mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    if (ncells_owned != ncolumns_local * ncells_per_column) {
      Errors::Message msg;
      msg << "InputOutputVerticallyStructuredNetCDF: this rank's owned cell count (" << ncells_owned
          << ") is not ncolumns_local (" << ncolumns_local << ") * ncells_per_column ("
          << ncells_per_column << ") -- columns do not partition "
          << "this rank's owned cells evenly.";
      Exceptions::amanzi_throw(msg);
    }
  }

  int ncells_global = mesh->getMap(AmanziMesh::Entity_kind::CELL, false)->getGlobalNumElements();
  if (ncells_global != ncolumns_global * ncells_per_column) {
    Errors::Message msg;
    msg << "InputOutputVerticallyStructuredNetCDF: mesh's total cell count (" << ncells_global
        << ") is not ncolumns_global (" << ncolumns_global << ") * this rank's ncells_per_column ("
        << ncells_per_column << ") -- columns are not uniform height mesh-wide.";
    Exceptions::amanzi_throw(msg);
  }
}

} // namespace InputOutputVerticallyStructuredNetCDF_Helpers


//
// InputVerticallyStructuredNetCDF
//
InputVerticallyStructuredNetCDF::InputVerticallyStructuredNetCDF(
  const Comm_ptr_type& comm,
  const std::string& filename,
  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  : file_(std::make_unique<FileNetCDF>(comm, filename, FILE_READONLY)), info_(mesh)
{}


namespace {

bool
hasSubfieldNames(const Teuchos::ParameterList& attrs, std::size_t count)
{
  return attrs.isParameter("subfieldnames") &&
         attrs.get<Teuchos::Array<std::string>>("subfieldnames").size() == count;
}

} // namespace


void
InputVerticallyStructuredNetCDF::read(const Teuchos::ParameterList& attrs,
                                      MultiVector_type& vec) const
{
  notImplemented_("read(MultiVector_type&)");
}


void
InputVerticallyStructuredNetCDF::readRange(
  const Teuchos::ParameterList& attrs,
  std::size_t time_start,
  std::size_t time_count,
  Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> buf) const
{
  int ndofs = buf.extent(1);
  AMANZI_ASSERT(buf.span_is_contiguous());

  FileNetCDF::TimeRange time_chunk{ time_start, time_start + time_count };

  if (info_.is_columnar) {
    AMANZI_ASSERT((int)buf.extent(2) == info_.ncolumns_local * info_.ncells_per_column);
    if (ndofs == 1) {
      // Reinterpret buf's memory as [time_count, ncolumns_local,
      // ncells_per_column]; splitting the fastest-varying axis in two
      // changes no element's address under LayoutRight.
      Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> columnar(
        buf.data(), time_count, info_.ncolumns_local, info_.ncells_per_column);
      file_->readDistributedView(
        attrs.name(), columnar, 1, { "time", "ncolumn", "ncell_per_column" }, time_chunk);
    } else {
      Kokkos::View<double****, Kokkos::LayoutRight, Kokkos::HostSpace> columnar(
        buf.data(), time_count, ndofs, info_.ncolumns_local, info_.ncells_per_column);
      file_->readDistributedView(attrs.name(),
                                 columnar,
                                 2,
                                 { "time", "ndofs_" + attrs.name(), "ncolumn", "ncell_per_column" },
                                 time_chunk);
    }
  } else {
    // a surface mesh -- no columns, read flat/unraveled instead (same
    // packed shape family writeFlat_ writes, named "ncell" since this
    // method only ever reads CELL data).
    if (ndofs == 1) {
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> flat(
        buf.data(), time_count, buf.extent(2));
      file_->readDistributedView(attrs.name(), flat, 1, { "time", "ncell" }, time_chunk);
    } else {
      file_->readDistributedView(
        attrs.name(), buf, 2, { "time", "ndofs_" + attrs.name(), "ncell" }, time_chunk);
    }
  }
}


//
// OutputVerticallyStructuredNetCDF
//
OutputVerticallyStructuredNetCDF::OutputVerticallyStructuredNetCDF(
  Teuchos::ParameterList& plist,
  const Comm_ptr_type& comm,
  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  : comm_(comm), info_(mesh)
{
  filename_ = plist.get<std::string>("file name base") + ".nc";
  file_ = std::make_unique<FileNetCDF>(comm_, filename_, FILE_CREATE);
  time_index_ = -1;
}


void
OutputVerticallyStructuredNetCDF::createTimestep(double time, int cycle)
{
  cycle_ = cycle;
  time_ = time;
  ++time_index_;
}


void
OutputVerticallyStructuredNetCDF::finalizeTimestep()
{
  AMANZI_ASSERT(time_index_ >= 0); // createTimestep() must precede finalizeTimestep()
  file_->writeAttribute("time", time_, (std::size_t)time_index_);
}


std::string
OutputVerticallyStructuredNetCDF::getFilename(int cycle) const
{
  return filename_;
}


void
OutputVerticallyStructuredNetCDF::write(const Teuchos::ParameterList& attrs, const int& val) const
{
  if (time_index_ < 0) {
    file_->writeAttribute(attrs.name(), val);
  } else {
    file_->writeAttribute(attrs.name(), val, (std::size_t)time_index_);
  }
}


void
OutputVerticallyStructuredNetCDF::write(const Teuchos::ParameterList& attrs,
                                        const double& val) const
{
  if (time_index_ < 0) {
    file_->writeAttribute(attrs.name(), val);
  } else {
    file_->writeAttribute(attrs.name(), val, (std::size_t)time_index_);
  }
}


void
OutputVerticallyStructuredNetCDF::write(const Teuchos::ParameterList& attrs,
                                        const std::string& val) const
{
  file_->writeAttribute<std::string>(attrs.name(), val);
}


void
OutputVerticallyStructuredNetCDF::write(const Teuchos::ParameterList& attrs,
                                        const Vector_type& vec) const
{
  write_(attrs.name(), vec, attrs.get<AmanziMesh::Entity_kind>("location"));
}


void
OutputVerticallyStructuredNetCDF::write(const Teuchos::ParameterList& attrs,
                                        const IntVector_type& vec) const
{
  Errors::Message msg("OutputVerticallyStructuredNetCDF::write(IntVector_type) is not supported "
                      "by FileNetCDF (double data only).");
  Exceptions::amanzi_throw(msg);
}


void
OutputVerticallyStructuredNetCDF::write_(const std::string& varname,
                                         const MultiVector_type& vec,
                                         AmanziMesh::Entity_kind location) const
{
  if (location == AmanziMesh::Entity_kind::CELL && info_.is_columnar) {
    writeRavelled_(varname, vec);
  } else {
    writeFlat_(varname, vec, location);
  }
}


void
OutputVerticallyStructuredNetCDF::writeFlat_(const std::string& varname,
                                             const MultiVector_type& vec,
                                             AmanziMesh::Entity_kind location) const
{
  AMANZI_ASSERT(time_index_ >= 0); // createTimestep() must precede any field write
  std::size_t time_index = (std::size_t)time_index_;
  int ndofs = vec.getNumVectors();
  int n = (int)vec.getMap()->getLocalNumElements();
  auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadOnly);

  // named by location, not varname, so every field at the same location
  // shares one dimension on disk -- "n" + singular entity name, matching
  // UGRID's convention for mesh-entity-count dimensions (e.g. nMesh2_face).
  std::string loc_name = "n" + AmanziMesh::to_string(location);

  if (ndofs == 1) {
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> buf(
      "OutputVerticallyStructuredNetCDF::writeFlat_", n);
    for (int i = 0; i != n; ++i) buf(i) = vecv(i, 0);
    file_->writeDistributedView(varname, buf, 0, { "time", loc_name }, time_index);
  } else {
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf(
      "OutputVerticallyStructuredNetCDF::writeFlat_", ndofs, n);
    for (int j = 0; j != ndofs; ++j) {
      for (int i = 0; i != n; ++i) buf(j, i) = vecv(i, j);
    }
    file_->writeDistributedView(
      varname, buf, 1, { "time", "ndofs_" + varname, loc_name }, time_index);
  }
}


void
OutputVerticallyStructuredNetCDF::writeRavelled_(const std::string& varname,
                                                 const MultiVector_type& vec) const
{
  AMANZI_ASSERT(time_index_ >= 0); // createTimestep() must precede any field write
  std::size_t time_index = (std::size_t)time_index_;
  int ndofs = vec.getNumVectors();
  auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadOnly);

  if (ndofs == 1) {
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf(
      "OutputVerticallyStructuredNetCDF::writeRavelled_",
      info_.ncolumns_local,
      info_.ncells_per_column);
    for (int col = 0; col != info_.ncolumns_local; ++col) {
      for (int layer = 0; layer != info_.ncells_per_column; ++layer) {
        int local_id = col * info_.ncells_per_column + layer;
        buf(col, layer) = vecv(local_id, 0);
      }
    }
    file_->writeDistributedView(
      varname, buf, 0, { "time", "ncolumn", "ncell_per_column" }, time_index);
  } else {
    Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> buf(
      "OutputVerticallyStructuredNetCDF::writeRavelled_",
      ndofs,
      info_.ncolumns_local,
      info_.ncells_per_column);
    for (int j = 0; j != ndofs; ++j) {
      for (int col = 0; col != info_.ncolumns_local; ++col) {
        for (int layer = 0; layer != info_.ncells_per_column; ++layer) {
          int local_id = col * info_.ncells_per_column + layer;
          buf(j, col, layer) = vecv(local_id, j);
        }
      }
    }
    file_->writeDistributedView(
      varname, buf, 1, { "time", "ndofs_" + varname, "ncolumn", "ncell_per_column" }, time_index);
  }
}


void
OutputVerticallyStructuredNetCDF::write(const Teuchos::ParameterList& attrs,
                                        const MultiVector_type& vec) const
{
  int ndofs = vec.getNumVectors();
  AmanziMesh::Entity_kind location = attrs.get<AmanziMesh::Entity_kind>("location");

  if (hasSubfieldNames(attrs, ndofs)) {
    std::vector<std::string> names = OutputUtils::getNames(attrs, ndofs);
    for (int j = 0; j != ndofs; ++j) { write_(names[j], *vec.getVector(j), location); }
  } else {
    write_(attrs.name(), vec, location);
  }
}

} // namespace Amanzi
