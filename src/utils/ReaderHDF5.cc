/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! ReaderHDF5: simple reader for serial reads of HDF5 files.

#include <string>

#include "Teuchos_Array.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "errors.hh"
#include "dbc.hh"
#include "ReaderHDF5.hh"

namespace Amanzi {

ReaderHDF5::ReaderHDF5(const std::string& filename)
  : filename_(filename), file_(-1)
{
  htri_t ierr = H5Fis_hdf5(filename.c_str());
  if (ierr > 0) {
    file_ = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  } else {
    Errors::Message msg;
    msg << "ReaderHDF5: error opening file \"" << filename << "\" with READ_ONLY access.";
    Exceptions::amanzi_throw(msg);
  }
}


ReaderHDF5::~ReaderHDF5()
{
  if (file_ >= 0) H5Fclose(file_);
}


bool
ReaderHDF5::hasVariableOrGroup(const std::string& varname) const
{
  return (H5Lexists(file_, varname.c_str(), H5P_DEFAULT) > 0);
}


void
ReaderHDF5::read(const std::string& lvarname, Teuchos::Array<double>& vec, int index) const
{
  std::string varname = lvarname;
  if (index >= 0) varname = varname + "/" + std::to_string(index);

  hid_t dataset = H5Dopen(file_, varname.c_str(), H5P_DEFAULT);
  if (dataset < 0) {
    Errors::Message msg;
    msg << "ReaderHDF5: requested variable \"" << varname << "\" is not found in file \""
        << filename_ << "\"";
    Exceptions::amanzi_throw(msg);
  }

  hid_t dataspace = H5Dget_space(dataset);
  hssize_t size = H5Sget_simple_extent_npoints(dataspace);
  vec.resize(size);
  herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec.data());
  if (status) {
    Errors::Message msg;
    msg << "ReaderHDF5: error reading variable \"" << varname << "\" in file \"" << filename_
        << "\"";
    Exceptions::amanzi_throw(msg);
  }
  H5Dclose(dataset);
}

void
ReaderHDF5::read(const std::string& lvarname, Teuchos::Array<int>& vec, int index) const
{
  std::string varname = lvarname;
  if (index >= 0) varname = varname + "/" + std::to_string(index);

  hid_t dataset = H5Dopen(file_, varname.c_str(), H5P_DEFAULT);
  if (dataset < 0) {
    Errors::Message msg;
    msg << "ReaderHDF5: requested variable \"" << varname << "\" is not found in file \""
        << filename_ << "\"";
    Exceptions::amanzi_throw(msg);
  }

  hid_t dataspace = H5Dget_space(dataset);
  hssize_t size = H5Sget_simple_extent_npoints(dataspace);
  vec.resize(size);
  herr_t status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec.data());
  if (status) {
    Errors::Message msg;
    msg << "ReaderHDF5: error reading variable \"" << varname << "\" in file \"" << filename_
        << "\"";
    Exceptions::amanzi_throw(msg);
  }
  H5Dclose(dataset);
}

void
ReaderHDF5::read(const std::string& lvarname,
                 Teuchos::SerialDenseMatrix<int, double>& mat,
                 int index) const
{
  std::string varname = lvarname;
  if (index >= 0) varname = varname + "/" + std::to_string(index);

  hid_t dataset = H5Dopen(file_, varname.c_str(), H5P_DEFAULT);
  if (dataset < 0) {
    Errors::Message msg;
    msg << "ReaderHDF5: requested variable \"" << varname << "\" is not found in file \""
        << filename_ << "\"";
    Exceptions::amanzi_throw(msg);
  }

  hid_t dataspace = H5Dget_space(dataset);
  hsize_t dims[2];
  int rank = H5Sget_simple_extent_dims(dataspace, dims, NULL);
  if (rank != 2) {
    Errors::Message message("ReaderHDF5: error, dataset dimension is not equal to 2 ");
    Exceptions::amanzi_throw(message);
  }

  // NOTE: HDF5 is row-major, SerialDenseMatrix is column-major
  mat.shapeUninitialized(dims[1], dims[0]);
  herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mat.values());
  if (status) {
    std::string message("ReaderHDF5: read error");
    Exceptions::amanzi_throw(message);
  }
  H5Dclose(dataset);
}

} // namespace Amanzi
