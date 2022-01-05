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
FileHDF5::WriteAttribute(const std::string& attr_name_,
                         const std::string& h5path_, scalar_type value)
{
  if (!HasGroup(h5path_)) {
    CreateGroup(h5path_);
    CloseGroup();
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
FileHDF5::ReadAttribute(const std::string& attr_name_,
                        const std::string& h5path_)
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
FileHDF5::WriteVector(const std::string& var_name,
                      const Vector_type_<scalar_type>& vec)
{
  IODetails::DangerousString full_h5path(var_name);

  // parallelIO only deals with global lengths that are ints
  // Really check this!
  static_assert(std::is_same<Tpetra::global_size_t, std::size_t>::value,
                "ASCEMIO does not support vectors with global length larger "
                "than a (signed) int");
  Tpetra::global_size_t globallength_raw = vec.getGlobalLength();

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


#if 0
template<typename scalar_type>
void
FileHDF5::WriteMultiVectorBlock(const std::string& var_name, const MultiVector_type_<scalar_type>& vec)
{
  // parallelIO only deals with global lengths that are ints
  // Really check this!
  static_assert(std::is_same<Tpetra::global_size_t, std::size_t>::value,
                "ASCEMIO does not support vectors with global length larger than a (signed) int");
  Tpetra::global_size_t globallength_raw = vec.getGlobalLength();

  if (globallength_raw > std::numeric_limits<int>::max()) {
    Errors::Message message("FileHDF5::ASCEMIO does not support vectors with global length larger than a (signed) int.");
    throw(message);
  }

  // must manage both layouts -- MultiVector could be either!
  typedef decltype(vec.getLocalViewHost()) host_view_type;

  int globaldims[2];
  int localdims[2];

  if (std::is_same<Kokkos::LayoutLeft, typename host_view_type::array_layout>::value) {
    // LayoutLeft, so entities are fastest-varying, and we must flip the dimensions
    globaldims[1] = static_cast<int>(globallength_raw);
    globaldims[0] = static_cast<int>(vec.getNumVectors());

    localdims[1] = static_cast<int>(vec.getLocalLength());
    localdims[0] = globaldims[0];
  } else {
    // LayoutRight, so vectors are fastest-varying, and dimensions are the same
    globaldims[0] = static_cast<int>(globallength_raw);
    globaldims[1] = static_cast<int>(vec.getNumVectors());

    localdims[0] = static_cast<int>(vec.getLocalLength());
    localdims[1] = globaldims[0];
  }

  // THIS SHOULD WORK BUT FAILS.  ASCEMIO has some wierd issues, but it appears
  // to me that it ensures that each IO node takes a block of ROWS (or the
  // slower-varying dimension), which is convenient if the decomposed dimension
  // is the slower-varying dimension (basically each IO node would get all data
  // from a set of processes, correct for PETSc), but is crazy bad if the
  // decomposed dimension is the faster-varying dimension (each IO node would
  // get a single row that touches all processes, this is Trilinos ordering).
  IODetails::DangerousString full_h5path(var_name);
  parallelIO_write_dataset((void*) vec.get1dView().get(),
                           IODetails::PIO_DatatypeMap<scalar_type>::type,
                           2, globaldims, localdims,
                           data_file_, full_h5path.c_str(),
                           &IOgroup_, NONUNIFORM_CONTIGUOUS_WRITE);
}
#endif

template <typename scalar_type>
void
FileHDF5::WriteMultiVector(const std::vector<std::string>& var_paths,
                           const MultiVector_type_<scalar_type>& vec)
{
  if (var_paths.size() != vec.getNumVectors()) {
    Errors::Message message(
      "FileHDF5::WriteMultiVector requested with invalid number of names.");
    throw(message);
  }

  // parallelIO only deals with global lengths that are ints
  // Really check this!
  static_assert(std::is_same<Tpetra::global_size_t, std::size_t>::value,
                "ASCEMIO does not support vectors with global length larger "
                "than a (signed) int");
  Tpetra::global_size_t globallength_raw = vec.getGlobalLength();

  if (globallength_raw > std::numeric_limits<int>::max()) {
    Errors::Message message("FileHDF5::ASCEMIO does not support vectors with "
                            "global length larger than a (signed) int.");
    throw(message);
  }
  int globaldims[] = { static_cast<int>(globallength_raw), 1 };

  std::size_t locallength = vec.getLocalLength();
  int localdims[] = { static_cast<int>(locallength), 1 };

  for (std::size_t i = 0; i != vec.getNumVectors(); ++i) {
    IODetails::DangerousString full_h5path(var_paths[i]);

    auto vec_data = vec.getData(i);
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
}


//
// NOTE: this is somewhat dangerous -- Views are a LOCAL concept, but this is a
// parallel write, so we have work to do.
//
template <typename view_type>
void
FileHDF5::WriteView(const std::string& var_path,
                    Tpetra::global_size_t global_length_raw,
                    const view_type vec)
{
  Kokkos::View<double**, Kokkos::LayoutRight, DefaultHost>
    check_layout_right = vec;
  int localdims[] = { static_cast<int>(vec.extent(0)),
                      static_cast<int>(vec.extent(1)) };

  IODetails::DangerousString full_h5path(var_path);
  // Fix for gcc 9.4.0
  std::size_t local_length = vec.extent(0);
  auto global_length(local_length);
  Teuchos::reduceAll(
    *comm_, Teuchos::REDUCE_SUM, 1, &local_length, &global_length);
  int globaldims[] = { static_cast<int>(global_length), localdims[1] };

  parallelIO_write_dataset((void*)&vec(0, 0),
                           PIO_DOUBLE,
                           2,
                           globaldims,
                           localdims,
                           data_file_,
                           full_h5path.c_str(),
                           &IOgroup_,
                           NONUNIFORM_CONTIGUOUS_WRITE);
}


// //
// // NOTE: this is somewhat dangerous -- Views are a LOCAL concept, but this is
// a
// // parallel read, so we have work to do.
// //
// template<typename view_type>
// void
// FileHDF5::ReadView(const std::string& var_name, view_type vec)
// {
//   IODetails::DangerousString full_h5path(var_name);

//   Kokkos::View<double**, Kokkos::LayoutRight, DefaultHost>
//   check_layout_right = vec; int localdims[] = {
//   static_cast<int>(vec.extent(0)), static_cast<int>(vec.extent(1)) };

//   int ndims;
//   parallelIO_get_dataset_ndims(&ndims, data_file_, full_h5path.c_str(),
//   &IOgroup_); if (ndims != 2) {
//     Errors::Message message("FileHDF5::ReadView only deals with 2D views.");
//     throw(message);
//   }

//   int globaldims[ndims];
//   parallelIO_get_dataset_dims(globaldims, data_file_, full_h5path.c_str(),
//   &IOgroup_);

//   // check dims add up
//   auto local_length = vec.extent(0);
//   auto global_length(local_length);
//   comm_->reduceAll(Teuchos::REDUCE_SUM, 1, &local_length, &global_length);
//   if (globaldims[0] != global_length) {
//     Errors::Message message("FileHDF5::ReadView called with incompatible
//     global dimension views."); throw(message);
//   }

//   for (int i=0; i!=ndims; ++i) {
//     if (vec.extent(i) != localdims[i]) {
//       Errors::Message message;
//       message << "FileHDF5::ReadView called with view of incorrect extent
//       (got " << localdims[i] << ", expected " << vec.extent(i) << " in
//       dimension " << i << ")."; throw(message);
//     }
//   }

//   localdims[0] = static_cast<int>(vec.getLocalLength());
//   localdims[1] = globaldims[1];

//   {
//     auto vecv = vec.get1dViewNonConst();
//     parallelIO_read_dataset((void *) vecv.get(),
//                             IODetails::PIO_DatatypeMap<scalar_type>::type,
//                             ndims, globaldims, localdims,
//                             data_file_, full_h5path.c_str(),
//                             &IOgroup_, NONUNIFORM_CONTIGUOUS_READ);
//   }
// }


template <typename scalar_type>
void
FileHDF5::ReadVector(const std::string& var_name,
                     Vector_type_<scalar_type>& vec)
{
  IODetails::DangerousString full_h5path(var_name);

  int ndims;
  parallelIO_get_dataset_ndims(
    &ndims, data_file_, full_h5path.c_str(), &IOgroup_);
  if (ndims < 0) {
    Errors::Message message(
      "FileHDF5::ReadVector data has negative dimension.");
    throw(message);
  }

  int globaldims[ndims], localdims[ndims];
  parallelIO_get_dataset_dims(
    globaldims, data_file_, full_h5path.c_str(), &IOgroup_);
  if (vec.getGlobalLength() != globaldims[0]) {
    Errors::Message message;
    message << "FileHDF5::ReadVector with incorrect length (got "
            << globaldims[0] << ", expected " << vec.getGlobalLength() << ").";
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
FileHDF5::ReadMultiVector(const std::vector<std::string>& var_paths,
                          MultiVector_type_<scalar_type>& vec)
{
  auto vecv = vec.get2dViewNonConst();

  if (var_paths.size() != vec.getNumVectors()) {
    Errors::Message message(
      "FileHDF5::ReadMultiVector requested with invalid number of names.");
    throw(message);
  }

  for (std::size_t i = 0; i != vec.getNumVectors(); ++i) {
    IODetails::DangerousString full_h5path(var_paths[i]);

    int ndims;
    parallelIO_get_dataset_ndims(
      &ndims, data_file_, full_h5path.c_str(), &IOgroup_);
    if (ndims < 0) {
      Errors::Message message(
        "FileHDF5::ReadVector data has negative dimension.");
      throw(message);
    }

    int globaldims[ndims], localdims[ndims];
    parallelIO_get_dataset_dims(
      globaldims, data_file_, full_h5path.c_str(), &IOgroup_);
    if (vec.getGlobalLength() != globaldims[0]) {
      Errors::Message message;
      message << "FileHDF5::ReadVector with incorrect length (got "
              << globaldims[0] << ", expected " << vec.getGlobalLength()
              << ").";
      throw(message);
    }

    localdims[0] = static_cast<int>(vec.getLocalLength());
    localdims[1] = globaldims[1];

    {
      parallelIO_read_dataset((void*)vecv[i].get(),
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
}


#if 0

template<typename scalar_type>
void
FileHDF5::ReadMultiVectorBlock(const std::string& var_name, MultiVector_type_<scalar_type>& vec)
{
  IODetails::DangerousString full_h5path(var_name);

  int ndims;
  parallelIO_get_dataset_ndims(&ndims, data_file_, full_h5path.c_str(), &IOgroup_);
  if (ndims < 0) {
    Errors::Message message("FileHDF5::ReadVector data has negative dimension.");
    throw(message);
  }

  int globaldims_file[ndims], globaldims[ndims], localdims[ndims];
  parallelIO_get_dataset_dims(globaldims_file, data_file_, full_h5path.c_str(), &IOgroup_);

  typedef decltype(vec.getLocalViewHost()) host_view_type;
  if (std::is_same<Kokkos::LayoutLeft, typename host_view_type::array_layout>::value) {
    // LayoutLeft, so entities are fastest-varying, and we must flip the dimensions
    globaldims[1] = static_cast<int>(vec.getGlobalLength());
    globaldims[0] = static_cast<int>(vec.getNumVectors());

    localdims[1] = static_cast<int>(vec.getLocalLength());
    localdims[0] = globaldims[0];
  } else {
    // LayoutRight, so vectors are fastest-varying, and dimensions are the same
    globaldims[0] = static_cast<int>(vec.getGlobalLength());
    globaldims[1] = static_cast<int>(vec.getNumVectors());

    localdims[0] = static_cast<int>(vec.getLocalLength());
    localdims[1] = globaldims[0];
  }

  if (globaldims[0] != globaldims_file[0] || globaldims[1] != globaldims_file[1]) {
    // ERROR!
    if (globaldims[0] == globaldims_file[1] && globaldims[1] == globaldims_file[0]) {
      Errors::Message message("FileHDF5::ReadMultiVectorBlock: Shape of MultiVector is transpose of shape of file.  Cannot use files generated on other machines reliably!");
      throw(message);
    } else {
      Errors::Message message;
      message << "FileHDF5::ReadMultiVectorBlock: Shape of MultiVector is different from shape of file, expected " << globaldims[0] << "," << globaldims[1]
              << " but got " << globaldims_file[0] << "," << globaldims_file[1] << ".";
      throw(message);
    }
  }

  {
    auto vecv = vec.get1dViewNonConst();
    parallelIO_read_dataset((void *) vecv.get(),
                            IODetails::PIO_DatatypeMap<scalar_type>::type,
                            ndims, globaldims, localdims,
                            data_file_, full_h5path.c_str(),
                            &IOgroup_, NONUNIFORM_CONTIGUOUS_READ);
  }
}

#endif


//
// STD::STRING
//
template <>
void
FileHDF5::WriteAttribute<std::string>(const std::string& attr_name,
                                      const std::string& h5path,
                                      std::string value);

template <>
std::string
FileHDF5::ReadAttribute<std::string>(const std::string& attr_name,
                                     const std::string& h5path);

} // namespace Amanzi
