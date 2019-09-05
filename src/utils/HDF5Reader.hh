/*
  Utils

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Interface for grabbing vectors of data from HDF5 files for use with i/o.
*/

#ifndef AMANZI_HDF5_READER_HH_
#define AMANZI_HDF5_READER_HH_

#include <vector>

#define H5Gcreate_vers 2
#define H5Dcreate_vers 2
#define H5Gopen_vers 2
#define H5Dopen_vers 2
#include "hdf5.h"
#include "Epetra_SerialDenseMatrix.h"

#include "errors.hh"

#include <Kokkos_Core.hpp>

namespace Amanzi {

struct HDF5Reader {
 public:
  HDF5Reader(std::string filename) :
      filename_(filename) {
    htri_t ierr = H5Fis_hdf5(filename.c_str());
    if (ierr > 0) {
      file_ = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    } else {
      std::string header("HDF5Reader: error, invalid filename ");
      Errors::Message message(header+filename);
      Exceptions::amanzi_throw(message);
    }
  }

  ~HDF5Reader() {
    H5Fclose(file_);
  }

  void
  ReadData(std::string varname, Kokkos::View<double*>& vec) {
    //    char *h5path = new char[varname.size()+1];
    //    strcpy(h5path,varname.c_str());

    hid_t dataset = H5Dopen(file_, varname.c_str(), H5P_DEFAULT);
    hid_t dataspace = H5Dget_space(dataset);
    hssize_t size = H5Sget_simple_extent_npoints(dataspace);
    //vec.resize(size);
    Kokkos::resize(vec,size); 
    herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE,  H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, vec.data());
  }

  void
  ReadMatData(std::string varname, Kokkos::View<double**> &mat) {
    //    char *h5path = new char[varname.size()+1];
    //    strcpy(h5path,varname.c_str());
    hid_t dataset = H5Dopen(file_, varname.c_str(), H5P_DEFAULT );
    hid_t dataspace = H5Dget_space(dataset);
    hsize_t dims[2];
    int rank = H5Sget_simple_extent_dims(dataspace, dims, NULL);
    if ( rank != 2 ) {
      Errors::Message message("HDF5Reader: error, dataset dimension is not equal to 2 ");
      Exceptions::amanzi_throw(message);
    }
    Kokkos::resize(mat,dims[0],dims[1]); 
    Epetra_SerialDenseMatrix m; 
    m.Shape(dims[1],dims[0]);

    // Not readable directly, need to change row major / column major 
    //herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE,  H5S_ALL, H5S_ALL,
    //                        H5P_DEFAULT, mat.data());

    herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE,  H5S_ALL, H5S_ALL,
                            H5P_DEFAULT, &m[0][0]);
    
    for(int i = 0 ; i < mat.extent(0) ; ++i){
      for(int j = 0 ; j < mat.extent(1); ++j){
        mat(i,j) = m[i][j]; 
      }
    }
  }

 protected:
  std::string filename_;
  hid_t file_;
};

} // namespace Amanzi

#endif
