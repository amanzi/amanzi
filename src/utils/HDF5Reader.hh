/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Interface for reading vectors of data from HDF5 files, in serial only.

#ifndef AMANZI_HDF5_READER_HH_
#define AMANZI_HDF5_READER_HH_

#include <vector>

#define H5Gcreate_vers 2
#define H5Dcreate_vers 2
#define H5Gopen_vers 2
#define H5Dopen_vers 2
#include "hdf5.h"

#include "Teuchos_Array.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "errors.hh"

#include <Kokkos_Core.hpp>

namespace Amanzi {

struct HDF5Reader {
 public:
  HDF5Reader(std::string filename) : filename_(filename)
  {
    htri_t ierr = H5Fis_hdf5(filename.c_str());
    if (ierr > 0) {
      file_ = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    } else {
      std::string header("HDF5Reader: error, invalid filename ");
      Errors::Message message(header + filename);
      Exceptions::amanzi_throw(message);
    }
  }

  ~HDF5Reader() { H5Fclose(file_); }

  void ReadData(std::string varname, Teuchos::Array<double>& arr)
  {
    hid_t dataset = H5Dopen(file_, varname.c_str(), H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);
    hssize_t size = H5Sget_simple_extent_npoints(dataspace);

    arr.resize(size);
    herr_t status = H5Dread(
      dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr.data());
  }
  void ReadData(std::string varname, std::vector<double>& arr)
  {
    hid_t dataset = H5Dopen(file_, varname.c_str(), H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);
    hssize_t size = H5Sget_simple_extent_npoints(dataspace);

    arr.resize(size);
    herr_t status = H5Dread(
      dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arr.data());
  }
  void ReadData(std::string varname, Kokkos::View<double*,Kokkos::HostSpace>& vec)
  {
    hid_t dataset = H5Dopen(file_, varname.c_str(), H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);
    hssize_t size = H5Sget_simple_extent_npoints(dataspace);
    Kokkos::resize(vec, size);
    herr_t status = H5Dread(
      dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vec.data());
  }

  void ReadMatData(std::string varname,
                   Teuchos::SerialDenseMatrix<std::size_t, double>& mat)
  {
    hid_t dataset = H5Dopen(file_, varname.c_str(), H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);
    hsize_t dims[2];
    int ndims = H5Sget_simple_extent_dims(dataspace, dims, NULL);
    if (ndims != 2) {
      Errors::Message message("HDF5Reader: dataset dimension is not 2.");
      Exceptions::amanzi_throw(message);
    }
    mat.shape(dims[1], dims[0]);
    herr_t status = H5Dread(
      dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, mat.values());
  }
  void ReadMatData(std::string varname, Kokkos::View<double**,Kokkos::HostSpace>& mat)
  {
    Teuchos::SerialDenseMatrix<std::size_t, double> m;
    ReadMatData(varname, m);

    Kokkos::resize(mat, m.numCols(), m.numRows());
    for (int i = 0; i < mat.extent(0); ++i) {
      for (int j = 0; j < mat.extent(1); ++j) { mat(i, j) = m[i][j]; }
    }
  }

 protected:
  std::string filename_;
  hid_t file_;
};

} // namespace Amanzi

#endif
