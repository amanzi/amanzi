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

#include "Tpetra_ConfigDefs.hpp"

#include "dbc.hh"
#include "errors.hh"
#include "AmanziVector.hh"

namespace Amanzi {

namespace IODetails {

// horrendous class to deal with const issues in ascemio
class DangerousString {
 public:
  DangerousString(std::string str) : str_(std::move(str)) {}
  DangerousString(const char* chars) : DangerousString(std::string(chars)) {}
  char* c_str() { return const_cast<char*>(str_.c_str()); }

 private:
  std::string str_;
};

// a compile-time map from type to ASCEMIO::datatype_t
template <typename T>
struct PIO_DatatypeMap;
template <>
struct PIO_DatatypeMap<int> {
  static const datatype_t type = PIO_INTEGER;
};
template <>
struct PIO_DatatypeMap<float> {
  static const datatype_t type = PIO_FLOAT;
};
template <>
struct PIO_DatatypeMap<double> {
  static const datatype_t type = PIO_DOUBLE;
};
template <>
struct PIO_DatatypeMap<long> {
  static const datatype_t type = PIO_LONG;
};
template <>
struct PIO_DatatypeMap<bool> {
  static const datatype_t type = PIO_BYTE;
};

} // namespace IODetails

template <typename scalar_type>
void
FileHDF5::writeAttribute(const std::string& attr_name_,
                         const std::string& h5path_,
                         scalar_type value)
{
  if (!hasGroup(h5path_)) {
    createGroup(h5path_);
    closeGroup();
  }
  IODetails::DangerousString attr_name(attr_name_);
  IODetails::DangerousString h5path(h5path_);
  parallelIO_write_simple_attr(attr_name.c_str(),
                               (void*)&value,
                               IODetails::PIO_DatatypeMap<scalar_type>::type,
                               data_file_,
                               h5path.c_str(),
                               &IOgroup_);
}


template <typename scalar_type>
scalar_type
FileHDF5::readAttribute(const std::string& attr_name_, const std::string& h5path_)
{
  IODetails::DangerousString attr_name(attr_name_);
  IODetails::DangerousString h5path(h5path_);

  scalar_type* loc_value;
  parallelIO_read_simple_attr(attr_name.c_str(),
                              reinterpret_cast<void**>(&loc_value),
                              IODetails::PIO_DatatypeMap<scalar_type>::type,
                              data_file_,
                              h5path.c_str(),
                              &IOgroup_);
  scalar_type value = *loc_value;
  free(loc_value);
  return value;
}


template <typename scalar_type>
void
FileHDF5::writeVector(const std::string& var_name, const Vector_type_<scalar_type>& vec)
{
  IODetails::DangerousString full_h5path(var_name);

  // parallelIO only deals with global lengths that are ints
  // Really check this!
  static_assert(std::is_same<Tpetra::global_size_t, std::size_t>::value,
                "ASCEMIO does not support vectors with global length larger "
                "than a (signed) int");
  GO globallength_raw = vec.getGlobalLength();

  if (globallength_raw > std::numeric_limits<int>::max()) {
    Errors::Message message("FileHDF5::ASCEMIO does not support vectors with "
                            "global length larger than a (signed) int.");
    throw(message);
  }
  int globaldims[] = { static_cast<int>(globallength_raw), 1 };

  std::size_t locallength = vec.getLocalLength();
  int localdims[] = { static_cast<int>(locallength), 1 };

  auto vec_data = vec.getData();
  parallelIO_write_dataset((void*)vec_data.get(),
                           IODetails::PIO_DatatypeMap<scalar_type>::type,
                           2,
                           globaldims,
                           localdims,
                           data_file_,
                           full_h5path.c_str(),
                           &IOgroup_,
                           NONUNIFORM_CONTIGUOUS_WRITE);
}


template <typename scalar_type>
void
FileHDF5::readVector(const std::string& var_name, Vector_type_<scalar_type>& vec)
{
  IODetails::DangerousString full_h5path(var_name);

  int ndims;
  parallelIO_get_dataset_ndims(&ndims, data_file_, full_h5path.c_str(), &IOgroup_);
  if (ndims < 0) {
    Errors::Message message("FileHDF5::readVector data has negative dimension.");
    throw(message);
  }

  int globaldims[ndims], localdims[ndims];
  parallelIO_get_dataset_dims(globaldims, data_file_, full_h5path.c_str(), &IOgroup_);
  if (vec.getGlobalLength() != globaldims[0]) {
    Errors::Message message;
    message << "FileHDF5::readVector with incorrect length (got " << globaldims[0] << ", expected "
            << vec.getGlobalLength() << ").";
    throw(message);
  }

  localdims[0] = static_cast<int>(vec.getLocalLength());
  localdims[1] = globaldims[1];

  {
    auto vecv = vec.get1dViewNonConst();
    parallelIO_read_dataset((void*)vecv.get(),
                            IODetails::PIO_DatatypeMap<scalar_type>::type,
                            ndims,
                            globaldims,
                            localdims,
                            data_file_,
                            full_h5path.c_str(),
                            &IOgroup_,
                            NONUNIFORM_CONTIGUOUS_READ);
  }
}


template <typename scalar_type>
void
FileHDF5::writeMultiVector(const std::string& var_name, const MultiVector_type_<scalar_type>& vec)
{
  writeView(var_name, vec.getLocalViewHost(Tpetra::Access::ReadOnly));
}


template <typename view_type>
void
FileHDF5::writeView(const std::string& var_path, view_type view)
{
  std::array<GO, 2> globaldims = { static_cast<GO>(view.extent(0)),
                                   static_cast<GO>(view.extent(1)) };
  Teuchos::reduceAll(*comm_, Teuchos::REDUCE_SUM, 2, globaldims.data(), globaldims.data());

  // ensure the shape is correct -- this is a PARALLEL write with decomposition in the 0th dimension
  AMANZI_ASSERT(globaldims[1] == static_cast<GO>(view.extent(1)) * comm_->getSize());
  globaldims[1] = view.extent(1);

  using const_scalar_type = typename view_type::const_value_type;
  using non_const_scalar_type = typename view_type::non_const_value_type;
  Kokkos::View<const_scalar_type**, Kokkos::LayoutRight, Kokkos::HostSpace> writeable_view;
  if constexpr (std::is_same<typename view_type::array_layout, Kokkos::LayoutRight>::value &&
                std::is_same<typename view_type::memory_space, Kokkos::HostSpace>::value) {
    writeable_view = view;
  } else {
    Kokkos::View<non_const_scalar_type**, Kokkos::LayoutRight, Kokkos::HostSpace> nc_writeable_view;
    Kokkos::resize(nc_writeable_view, view.extent(0), view.extent(1));
    Kokkos::deep_copy(nc_writeable_view, view);
    writeable_view = nc_writeable_view;
  }

  writeView(var_path, writeable_view, globaldims);
}


template <typename scalar_type>
void
FileHDF5::writeView(const std::string& var_path,
                    const Kokkos::View<scalar_type**, Kokkos::LayoutRight, Kokkos::HostSpace>& view,
                    std::array<GO, 2> dims)
{
  int globaldims[2], localdims[2];
  AMANZI_ASSERT(dims[0] < std::numeric_limits<int>::max());
  globaldims[0] = static_cast<int>(dims[0]);
  globaldims[1] = static_cast<int>(dims[1]);
  localdims[0] = static_cast<int>(view.extent(0));
  localdims[1] = static_cast<int>(view.extent(1));
  AMANZI_ASSERT(globaldims[1] == localdims[1]);

  IODetails::DangerousString full_h5path(var_path);

  using value_type = typename Kokkos::View<scalar_type**, Kokkos::LayoutRight, Kokkos::HostSpace>::
    non_const_value_type;
  parallelIO_write_dataset((void*)view.data(),
                           IODetails::PIO_DatatypeMap<value_type>::type,
                           2,
                           globaldims,
                           localdims,
                           data_file_,
                           full_h5path.c_str(),
                           &IOgroup_,
                           NONUNIFORM_CONTIGUOUS_WRITE);
}


template <typename scalar_type>
void
FileHDF5::readMultiVector(const std::string& var_name, MultiVector_type_<scalar_type>& vec)
{
  readView(var_name, vec.getLocalViewHost(Tpetra::Access::ReadWrite));
}


template <typename view_type>
void
FileHDF5::readView(const std::string& var_path, view_type view)
{
  std::array<GO, 2> globaldims = { static_cast<GO>(view.extent(0)),
                                   static_cast<GO>(view.extent(1)) };
  Teuchos::reduceAll(*comm_, Teuchos::REDUCE_SUM, 2, globaldims.data(), globaldims.data());

  // ensure the shape is correct -- this is a PARALLEL read with decomposition in the 1st dimension
  AMANZI_ASSERT(globaldims[1] == static_cast<GO>(view.extent(1)) * comm_->getSize());
  globaldims[1] = view.extent(1);

  using scalar_type = typename view_type::non_const_value_type;
  Kokkos::View<scalar_type**, Kokkos::LayoutRight, Kokkos::HostSpace> readable_view;
  bool deep_copy = false;
  if constexpr (std::is_same<typename view_type::array_layout, Kokkos::LayoutRight>::value &&
                std::is_same<typename view_type::memory_space, Kokkos::HostSpace>::value) {
    readable_view = view;
  } else {
    deep_copy = true;
    Kokkos::resize(readable_view, view.extent(0), view.extent(1));
  }
  readView(var_path, readable_view, globaldims);

  if (deep_copy) Kokkos::deep_copy(view, readable_view);
}


template <typename scalar_type>
void
FileHDF5::readView(const std::string& var_path,
                   Kokkos::View<scalar_type**, Kokkos::LayoutRight, Kokkos::HostSpace>& view,
                   std::array<GO, 2> dims)
{
  AMANZI_ASSERT(dims[0] < std::numeric_limits<int>::max());
  int globaldims[] = { static_cast<int>(dims[0]), static_cast<int>(dims[1]) };
  int localdims[] = { static_cast<int>(view.extent(0)), static_cast<int>(view.extent(1)) };
  AMANZI_ASSERT(globaldims[1] == localdims[1]);

  IODetails::DangerousString full_h5path(var_path);
  using value_type = typename Kokkos::View<scalar_type**, Kokkos::LayoutRight, Kokkos::HostSpace>::
    non_const_value_type;
  parallelIO_read_dataset((void*)view.data(),
                          IODetails::PIO_DatatypeMap<value_type>::type,
                          2,
                          globaldims,
                          localdims,
                          data_file_,
                          full_h5path.c_str(),
                          &IOgroup_,
                          NONUNIFORM_CONTIGUOUS_READ);
}


//
// STD::STRING
//
template <>
void
FileHDF5::writeAttribute<std::string>(const std::string& attr_name,
                                      const std::string& h5path,
                                      std::string value);

template <>
std::string
FileHDF5::readAttribute<std::string>(const std::string& attr_name, const std::string& h5path);

} // namespace Amanzi
