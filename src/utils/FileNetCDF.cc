/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Reads and writes NetCDF files for extruded-column (MPAS-style) meshes.

#include "dbc.hh"
#include "FileNetCDF.hh"

namespace Amanzi {

FileNetCDF::FileNetCDF(const Comm_ptr_type& comm, const std::string& filename, file_mode_t mode)
  : comm_(comm), filename_(filename), mode_(mode), ncid_(-1), file_open_(false)
{
  auto mpi_comm = Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm_);
  AMANZI_ASSERT(mpi_comm.get());
  mpi_comm_ = *mpi_comm->getRawMpiComm();
  openFile(mode);
}


FileNetCDF::~FileNetCDF()
{
  closeFile();
}


void
FileNetCDF::checkErr_(int ierr, const std::string& what) const
{
  if (ierr != NC_NOERR) {
    Errors::Message msg;
    msg << "FileNetCDF: error in " << what << " on file \"" << filename_
        << "\": " << nc_strerror(ierr);
    Exceptions::amanzi_throw(msg);
  }
}


void
FileNetCDF::openFile(file_mode_t mode)
{
  if (ncid_ >= 0) return;

  int ierr;
  if (mode == FILE_CREATE) {
    ierr =
      nc_create_par(filename_.c_str(), NC_NETCDF4 | NC_CLOBBER, mpi_comm_, MPI_INFO_NULL, &ncid_);
  } else {
    int nc_mode = (mode == FILE_READONLY) ? NC_NOWRITE : NC_WRITE;
    ierr = nc_open_par(filename_.c_str(), nc_mode | NC_NETCDF4, mpi_comm_, MPI_INFO_NULL, &ncid_);
  }
  checkErr_(ierr, "openFile(\"" + filename_ + "\")");
  file_open_ = true;
  dimids_.clear();
  varids_.clear();
}


void
FileNetCDF::closeFile()
{
  if (ncid_ >= 0) {
    int ierr = nc_close(ncid_);
    AMANZI_ASSERT(ierr == NC_NOERR);
    ncid_ = -1;
    file_open_ = false;
  }
}


//
// STD::STRING attribute specializations, mirroring FileHDF5's.
//
template <>
void
FileNetCDF::writeAttribute<std::string>(const std::string& attr_name, std::string value)
{
  int ierr = nc_put_att_text(ncid_, NC_GLOBAL, attr_name.c_str(), value.size(), value.c_str());
  checkErr_(ierr, "writeAttribute<std::string>(\"" + attr_name + "\")");
}


template <>
std::string
FileNetCDF::readAttribute<std::string>(const std::string& attr_name)
{
  nc_type xtype;
  size_t len;
  int ierr = nc_inq_att(ncid_, NC_GLOBAL, attr_name.c_str(), &xtype, &len);
  checkErr_(ierr, "readAttribute<std::string> inq(\"" + attr_name + "\")");

  std::string value(len, '\0');
  ierr = nc_get_att_text(ncid_, NC_GLOBAL, attr_name.c_str(), &value[0]);
  checkErr_(ierr, "readAttribute<std::string>(\"" + attr_name + "\")");
  return value;
}


int
FileNetCDF::getOrCreateDim_(const std::string& dimname, size_t len, bool unlimited)
{
  auto it = dimids_.find(dimname);
  if (it != dimids_.end()) return it->second;

  int dimid;
  int ierr = nc_inq_dimid(ncid_, dimname.c_str(), &dimid);
  if (ierr == NC_NOERR) {
    dimids_[dimname] = dimid;
    return dimid;
  }

  // not found -- define it (only valid if the file is writable and still in
  // define mode, which NetCDF-4 keeps open implicitly).
  ierr = nc_def_dim(ncid_, dimname.c_str(), unlimited ? NC_UNLIMITED : len, &dimid);
  checkErr_(ierr, "getOrCreateDim_(\"" + dimname + "\")");
  dimids_[dimname] = dimid;
  return dimid;
}


//
// Ensures a variable with the given on-disk dimension list exists.  dims is
// outermost-first, each entry's length already the dimension's full global
// length (computed by the caller -- see NetCDFDetails::buildDimsAndSlab,
// FileNetCDF_impl.hh).  dimids_out is populated with the resulting
// dimension ids, same order as dims.
//
int
FileNetCDF::getOrCreateVariable_(const std::string& varname,
                                 const std::vector<DimSpec>& dims,
                                 std::vector<int>& dimids_out,
                                 nc_type type)
{
  dimids_out.clear();
  for (const auto& dim : dims)
    dimids_out.push_back(getOrCreateDim_(dim.name, dim.length, dim.unlimited));

  auto it = varids_.find(varname);
  if (it != varids_.end()) return it->second;

  int varid;
  int ierr = nc_inq_varid(ncid_, varname.c_str(), &varid);
  if (ierr != NC_NOERR) {
    ierr =
      nc_def_var(ncid_, varname.c_str(), type, (int)dimids_out.size(), dimids_out.data(), &varid);
    checkErr_(ierr, "getOrCreateVariable_(\"" + varname + "\")");
  }

  // parallel access: collective (every rank participates in every
  // collective I/O call, even if its local piece is empty).
  // NC_VAR_PAR_ACCESS is a property of this open file handle (ncid_), not
  // persisted on disk -- so this must be set every time the variable is
  // first seen by THIS handle, not only when it is newly defined (varids_
  // is per-instance, so this runs at most once per (instance, varname), on
  // reopen included).
  ierr = nc_var_par_access(ncid_, varid, NC_COLLECTIVE);
  checkErr_(ierr, "nc_var_par_access(\"" + varname + "\")");

  varids_[varname] = varid;
  return varid;
}

} // namespace Amanzi
