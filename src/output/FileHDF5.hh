//! Reads and writes HDF5 files.

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

/*
  Reads and writes HDF5 files via parallelIO library.
*/

#ifndef AMANZI_FILE_HDF5_HH_
#define AMANZI_FILE_HDF5_HH_


extern "C" {
#include "hdf5.h"
#include "parallelIO.h"
};

#include "Tpetra_ConfigDefs.hpp"
#include "AmanziTypes.hh"

namespace Amanzi {

namespace AmanziMesh {
class Mesh;
}


class FileHDF5 {
 public:
  FileHDF5(const Comm_ptr_type& comm, const std::string& filename, file_mode_t mode);
  ~FileHDF5();

  void OpenFile(file_mode_t mode=FILE_READWRITE);
  void CloseFile();
  void CreateGroup(const std::string& h5path);
  

  template<typename scalar_type>
  void WriteAttribute(scalar_type value,
                 const std::string& attr_name, const std::string& h5path="/");

  template<typename scalar_type>
  void WriteAttributeArray(scalar_type const * const values, std::size_t count,
                      const std::string& attr_name, const std::string& h5path="/");

  template<typename scalar_type>
  scalar_type ReadAttribute(const std::string& attr_name, const std::string& h5path="/");

  template<typename scalar_type>
  std::vector<scalar_type> ReadAttributeArray(std::size_t count, const std::string& attr_name,
          const std::string& h5path="/");

  void WriteVector(const Map_type& vec, const std::string& var_path);
  template<typename scalar_type>
  void WriteVector(const Vector_type_<scalar_type>& vec, const std::string& var_path);
  template<typename view_type>
  void WriteView(const view_type vec, Tpetra::global_size_t global_length,
                 const std::string& var_path);
  template<typename scalar_type>
  void WriteMultiVectorBlock(const MultiVector_type_<scalar_type>& vec,
                        const std::string& var_name);
  template<typename scalar_type>
  void WriteMultiVector(const MultiVector_type_<scalar_type>& vec,
                        const std::vector<std::string>& var_paths);

  template<typename scalar_type>
  void ReadVector(Vector_type_<scalar_type>& vec, const std::string& var_name);
  
  template<typename scalar_type>
  void ReadMultiVector(MultiVector_type_<scalar_type>& vec,
                       const std::vector<std::string>& var_name);
  template<typename scalar_type>
  void ReadMultiVectorBlock(MultiVector_type_<scalar_type>& vec,
                            const std::string& var_name);
  

 private:
  iogroup_conf_t IOconfig_;
  iogroup_t IOgroup_;
  hid_t data_file_;
  std::vector<int> open_groups_;

  std::string filename_;
};  


} // namespace Amanzi


#include "FileHDF5_impl.hh"


#endif
