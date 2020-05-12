/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Reads and writes HDF5 files.

/*
  Reads and writes HDF5 files via parallelIO library.
*/

#ifndef AMANZI_FILE_HDF5_HH_
#define AMANZI_FILE_HDF5_HH_


extern "C"
{
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
  FileHDF5(const Comm_ptr_type& comm, const std::string& filename,
           file_mode_t mode);
  ~FileHDF5();

  void OpenFile(file_mode_t mode = FILE_READWRITE);
  void CloseFile();
  void CreateGroup(const std::string& h5path);
  void CloseGroup();
  bool HasGroup(const std::string& h5path) const;


  template <typename scalar_type>
  void WriteAttribute(const std::string& attr_name,
                      const std::string& attr_path, scalar_type value);

  template <typename scalar_type>
  scalar_type
  ReadAttribute(const std::string& attr_name, const std::string& h5path = "/");

  // write Maps
  void WriteVector(const std::string& var_path, const Map_type& vec);

  // read/write Vector_type_<>
  template <typename scalar_type>
  void WriteVector(const std::string& var_path,
                   const Vector_type_<scalar_type>& vec);

  template <typename scalar_type>
  void ReadVector(const std::string& var_name, Vector_type_<scalar_type>& vec);

  // write Kokkos::View<double**> objects.  Note this looks more general
  // because it takes a more general template argument, but really it is 2D
  // doubles only!
  template <typename view_type>
  void WriteView(const std::string& var_name,
                 Tpetra::global_size_t global_length, const view_type vec);

  // write MultiVector_type_<> in a single block (fastest)
  template <typename scalar_type>
  void WriteMultiVectorBlock(const std::string& var_name,
                             const MultiVector_type_<scalar_type>& vec);

  template <typename scalar_type>
  void ReadMultiVector(const std::vector<std::string>& var_name,
                       MultiVector_type_<scalar_type>& vec);

  // write MultiVector_type_<> by vector (currently how vis is done, fix me!)
  template <typename scalar_type>
  void WriteMultiVector(const std::vector<std::string>& var_paths,
                        const MultiVector_type_<scalar_type>& vec);

  template <typename scalar_type>
  void ReadMultiVectorBlock(const std::string& var_name,
                            MultiVector_type_<scalar_type>& vec);


 private:
  Comm_ptr_type comm_;

  iogroup_conf_t IOconfig_;
  iogroup_t IOgroup_;
  hid_t data_file_;
  std::vector<int> open_groups_;
  std::vector<std::string> all_groups_;

  std::string filename_;
};


} // namespace Amanzi


#include "FileHDF5_impl.hh"


#endif
