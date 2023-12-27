/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! reads and writes HDF5 files.
/*
  reads and writes HDF5 files via parallelIO library.
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

class FileHDF5 {
 public:
  FileHDF5(const Comm_ptr_type& comm, const std::string& filename, file_mode_t mode);
  ~FileHDF5();

  void openFile(file_mode_t mode = FILE_READWRITE);
  void closeFile();
  void createGroup(const std::string& h5path);
  void closeGroup();
  bool hasGroup(const std::string& h5path) const;

  template <typename scalar_type>
  void
  writeAttribute(const std::string& attr_name, const std::string& attr_path, scalar_type value);

  template <typename scalar_type>
  scalar_type readAttribute(const std::string& attr_name, const std::string& h5path = "/");

  // read/write Vector_type_<>
  template <typename scalar_type>
  void writeVector(const std::string& var_path, const Vector_type_<scalar_type>& vec);

  template <typename scalar_type>
  void readVector(const std::string& var_name, Vector_type_<scalar_type>& vec);

  template <typename scalar_type>
  void writeMultiVector(const std::string& var_path, const MultiVector_type_<scalar_type>& vec);

  template <typename scalar_type>
  void readMultiVector(const std::string& var_name, MultiVector_type_<scalar_type>& vec);

  template <typename view_type>
  void writeView(const std::string& var_path, view_type view);

  template <typename scalar_type>
  void writeView(const std::string& var_path,
                 const Kokkos::View<scalar_type**, Kokkos::LayoutRight, Kokkos::HostSpace>& view,
                 std::array<GO, 2> global_dims);

  template <typename view_type>
  void readView(const std::string& var_path, view_type vec);

  template <typename scalar_type>
  void readView(const std::string& var_path,
                Kokkos::View<scalar_type**, Kokkos::LayoutRight, Kokkos::HostSpace>& vec,
                std::array<GO, 2> global_dims);

  template <typename view_type>
  void readView(const std::string& var_path, const view_type vec, std::array<int, 2> global_dims);

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
