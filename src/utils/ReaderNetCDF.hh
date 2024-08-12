/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! NetCDFReader: simple reader for serial reads of netCDF4 files.
#pragma once

#include <string>

extern "C"
{
#include "netcdf.h"
};

#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "Key.hh"
#include "errors.hh"
#include "dbc.hh"
#include "Reader.hh"

namespace Amanzi {

class NetCDFReader : public Reader {
 public:
  explicit NetCDFReader(const std::string& filename);
  ~NetCDFReader();

  bool hasVariableOrGroup(const std::string& varname) const override;

  virtual void read(const std::string& varname, Teuchos::SerialDenseVector<int, double>& vec, int index = -1) const override {
    read_(varname, vec, index);
  }
  virtual void read(const std::string& varname, Teuchos::SerialDenseVector<int, int>& vec, int index = -1) const override {
    read_(varname, vec, index);
  }

  virtual void read(const std::string& varname, Teuchos::SerialDenseMatrix<int, double>& mat, int index = -1) const override {
    read_(varname, mat, index);
  }

  virtual void read(const std::string& varname, Teuchos::SerialDenseMatrix<int, int>& mat, int index = -1) const override {
    read_(varname, mat, index);
  }

 protected:
  // NOTE: intentional copy made here
  std::pair<int, int> findVarOrGroup_(std::string varname) const;

  template<typename Scalar>
  void read_(const std::string& varname, Teuchos::SerialDenseVector<int, Scalar>& vec, int index = -1) const;

  template<typename Scalar>
  void read_(const std::string& varname, Teuchos::SerialDenseMatrix<int, Scalar>& vec, int index = -1) const;

  std::string filename_;
  int file_;
};


template<typename Scalar>
void
NetCDFReader::read_(const std::string& varname, Teuchos::SerialDenseVector<int, Scalar>& vec, int index) const
{
  std::cout << "Reading " << varname << " index " << index << std::endl;
  auto [ncid, varid] = findVarOrGroup_(varname);
  if (varid < 0) {
    Errors::Message msg;
    msg << "NetCDFReader: error opening variable \"" << varname << "\" in file \"" << filename_ << "\"";
    Exceptions::amanzi_throw(msg);
  }

  // figure out the shape of var
  int ndims;
  int dim_ids[NC_MAX_VAR_DIMS];
  size_t dims[NC_MAX_VAR_DIMS];
  int ierr = nc_inq_var(ncid, varid, nullptr, nullptr, &ndims, dim_ids, nullptr);
  AMANZI_ASSERT(!ierr);
  for (int i = 0; i != ndims; ++i) {
    ierr = nc_inq_dim(ncid, dim_ids[i], nullptr, &dims[i]);
    AMANZI_ASSERT(!ierr);
  }

  if (index < 0) {
    if (ndims != 1) {
      Errors::Message msg;
      msg << "NetCDFReader: incorrect number of dimensions, " << ndims << ", for variable \"" << varname
          << "\" in file \"" << filename_ << "\" -- expected 1";
      Exceptions::amanzi_throw(msg);
    }
    vec.sizeUninitialized(dims[0]);
    ierr = nc_get_var(ncid, varid, vec.values());
  } else {
    if (ndims != 2) {
      Errors::Message msg;
      msg << "NetCDFReader: incorrect number of dimensions, " << ndims << ", for variable \"" << varname
          << "\" in file \"" << filename_ << "\" -- expected 2";
      Exceptions::amanzi_throw(msg);
    }
    if (index >= dims[0]) {
      Errors::Message msg;
      msg << "NetCDFReader: cannot read index " << index << " for variable \"" << varname
          << "\" in file \"" << filename_ << "\" as dimension is of length " << dims[0];
      Exceptions::amanzi_throw(msg);
    }
    vec.sizeUninitialized(dims[1]);

    size_t start[2] = { (size_t) index, 0 };
    size_t count[2] = { (size_t) 1, dims[1] };
    ierr = nc_get_vara(ncid, varid, start, count, vec.values());
    AMANZI_ASSERT(!ierr);
  }
}


template<typename Scalar>
void
NetCDFReader::read_(const std::string& varname, Teuchos::SerialDenseMatrix<int, Scalar>& mat, int index) const
{
  std::cout << "Reading Mat " << varname << " index " << index << std::endl;

  auto [ncid, varid] = findVarOrGroup_(varname);
  if (varid < 0) {
    Errors::Message msg;
    msg << "NetCDFReader: error opening variable \"" << varname << "\" in file \"" << filename_ << "\"";
    Exceptions::amanzi_throw(msg);
  }

  // figure out the shape of var
  int ndims;
  int dim_ids[NC_MAX_VAR_DIMS];
  size_t dims[NC_MAX_VAR_DIMS];
  int ierr = nc_inq_var(ncid, varid, nullptr, nullptr, &ndims, dim_ids, nullptr);
  AMANZI_ASSERT(!ierr);
  for (int i = 0; i != ndims; ++i) {
    ierr = nc_inq_dim(ncid, dim_ids[i], nullptr, &dims[i]);
    AMANZI_ASSERT(!ierr);
  }

  if (index < 0) {
    if (ndims != 2) {
      Errors::Message msg;
      msg << "NetCDFReader: incorrect number of dimensions, " << ndims << ", for variable \"" << varname
          << "\" in file \"" << filename_ << "\" -- expected 2";
      Exceptions::amanzi_throw(msg);
    }

    // NOTE: HDF5 is row-major, SerialDenseMatrix is column-major
    mat.shapeUninitialized(dims[1], dims[0]);
    ierr = nc_get_var(ncid, varid, mat.values());

  } else {
    if (ndims != 3) {
      Errors::Message msg;
      msg << "NetCDFReader: incorrect number of dimensions, " << ndims << ", for variable \"" << varname
          << "\" in file \"" << filename_ << "\" -- expected 3";
      Exceptions::amanzi_throw(msg);
    }
    if (index >= dims[0]) {
      Errors::Message msg;
      msg << "NetCDFReader: cannot read index " << index << " for variable \"" << varname
          << "\" in file \"" << filename_ << "\" as dimension is of length " << dims[0];
      Exceptions::amanzi_throw(msg);
    }

    // NOTE: HDF5 is row-major, SerialDenseMatrix is column-major
    mat.shapeUninitialized(dims[2], dims[1]);

    size_t start[3] = { (size_t) index, 0, 0 };
    size_t count[3] = { (size_t) 1, dims[1], dims[2] };
    ierr = nc_get_vara(ncid, varid, start, count, mat.values());
    AMANZI_ASSERT(!ierr);
  }
}

} // namespace Amanzi
