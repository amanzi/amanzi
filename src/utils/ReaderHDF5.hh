/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! HDF5Reader: simple reader for serial reads of HDF5 files.
#pragma once

#include <string>

extern "C"
{
#define H5Gcreate_vers 2
#define H5Dcreate_vers 2
#define H5Gopen_vers 2
#define H5Dopen_vers 2
#include "hdf5.h"
};

#include "dbc.hh"
#include "Reader.hh"

namespace Amanzi {

class HDF5Reader : public Reader {
 public:
  explicit HDF5Reader(const std::string& filename);
  ~HDF5Reader();

  bool hasVariableOrGroup(const std::string& varname) const override;

  void read(const std::string& varname, Teuchos::SerialDenseVector<int, double>& vec, int index = -1) const override;
  void read(const std::string& varname, Teuchos::SerialDenseVector<int, int>& vec, int index = -1) const override;
  void read(const std::string& varname, Teuchos::SerialDenseMatrix<int, double>& vec, int index = -1) const override;
  void read(const std::string& varname, Teuchos::SerialDenseMatrix<int, int>& vec, int index = -1) const override {
    AMANZI_ASSERT(false); // not implemented
  }

 protected:
  std::string filename_;
  hid_t file_;
};

} // namespace Amanzi
