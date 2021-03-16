/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

*/
//! HDF5Reader: simple reader for serial reads of HDF5 files.

#pragma once

#include <vector>

#define H5Gcreate_vers 2
#define H5Dcreate_vers 2
#define H5Gopen_vers 2
#define H5Dopen_vers 2
#include "hdf5.h"
#include "Epetra_SerialDenseMatrix.h"

#include "errors.hh"

namespace Amanzi {

class HDF5Reader {
 public:
  explicit HDF5Reader(const std::string& filename);
  ~HDF5Reader();

  void ReadData(std::string varname, std::vector<double>& vec);
  void ReadData(std::string varname, std::vector<int>& vec);
  void ReadMatData(std::string varname, Epetra_SerialDenseMatrix &mat);

 protected:
  std::string filename_;
  hid_t file_;

};

} // namespace Amanzi
