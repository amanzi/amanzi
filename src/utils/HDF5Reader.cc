/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! HDF5Reader: simple reader for serial reads of HDF5 files.
#include "HDF5Reader.hh"

namespace Amanzi {

HDF5Reader::HDF5Reader(const std::string& filename) : filename_(filename), file_(-1)
{
  htri_t ierr = H5Fis_hdf5(filename.c_str());
  if (ierr > 0) {
    file_ = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  } else {
    Errors::Message msg;
    msg << "HDF5Reader: error opening file \"" << filename << "\" with READ_ONLY access.";
    Exceptions::amanzi_throw(msg);
  }
}

HDF5Reader::~HDF5Reader()
{
  if (file_ >= 0) H5Fclose(file_);
}

bool
HDF5Reader::CheckVariableName(const std::string& varname)
{
  return (H5Lexists(file_, varname.c_str(), H5P_DEFAULT) > 0);
}

void
HDF5Reader::ReadData(std::string varname, std::vector<double>& vec)
{
  hid_t dataset = H5Dopen(file_, varname.c_str(), H5P_DEFAULT);
  if (dataset < 0) {
    Errors::Message msg;
    msg << "HDF5Reader: requested variable \"" << varname << "\" is not found in file \""
        << filename_ << "\"";
    Exceptions::amanzi_throw(msg);
  }

  hid_t dataspace = H5Dget_space(dataset);
  hssize_t size = H5Sget_simple_extent_npoints(dataspace);
  vec.resize(size);
  herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec.data());
  if (status) {
    Errors::Message msg;
    msg << "HDF5Reader: error reading variable \"" << varname << "\" in file \"" << filename_
        << "\"";
    Exceptions::amanzi_throw(msg);
  }
  H5Dclose(dataset);
}

void
HDF5Reader::ReadData(std::string varname, std::vector<int>& vec)
{
  hid_t dataset = H5Dopen(file_, varname.c_str(), H5P_DEFAULT);
  if (dataset < 0) {
    Errors::Message msg;
    msg << "HDF5Reader: requested variable \"" << varname << "\" is not found in file \""
        << filename_ << "\"";
    Exceptions::amanzi_throw(msg);
  }

  hid_t dataspace = H5Dget_space(dataset);
  hssize_t size = H5Sget_simple_extent_npoints(dataspace);
  vec.resize(size);
  herr_t status = H5Dread(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec.data());
  if (status) {
    Errors::Message msg;
    msg << "HDF5Reader: error reading variable \"" << varname << "\" in file \"" << filename_
        << "\"";
    Exceptions::amanzi_throw(msg);
  }
  H5Dclose(dataset);
}

void
HDF5Reader::ReadMatData(std::string varname, Epetra_SerialDenseMatrix& mat)
{
  hid_t dataset = H5Dopen(file_, varname.c_str(), H5P_DEFAULT);
  if (dataset < 0) {
    Errors::Message msg;
    msg << "HDF5Reader: requested variable \"" << varname << "\" is not found in file \""
        << filename_ << "\"";
    Exceptions::amanzi_throw(msg);
  }


  hid_t dataspace = H5Dget_space(dataset);
  hsize_t dims[2];
  int rank = H5Sget_simple_extent_dims(dataspace, dims, NULL);
  if (rank != 2) {
    Errors::Message message("HDF5Reader: error, dataset dimension is not equal to 2 ");
    Exceptions::amanzi_throw(message);
  }
  mat.Shape(dims[1], dims[0]);
  herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &mat[0][0]);
  if (status) {
    std::string message("HDF5Reader: read error");
    Exceptions::amanzi_throw(message);
  }
  H5Dclose(dataset);
}

} // namespace Amanzi
