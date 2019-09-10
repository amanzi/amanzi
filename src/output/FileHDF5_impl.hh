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
template<typename T> struct PIO_DatatypeMap;
template<> struct PIO_DatatypeMap<int> { static const datatype_t type = PIO_INTEGER; };
template<> struct PIO_DatatypeMap<float> { static const datatype_t type = PIO_FLOAT; };
template<> struct PIO_DatatypeMap<double> { static const datatype_t type = PIO_DOUBLE; };
template<> struct PIO_DatatypeMap<long> { static const datatype_t type = PIO_LONG; };
template<> struct PIO_DatatypeMap<bool> { static const datatype_t type = PIO_BYTE; };

} // namespace

template<typename scalar_type>
void
FileHDF5::WriteAttribute(scalar_type value,
               const std::string& attr_name_, const std::string& h5path_)
{
  IODetails::DangerousString attr_name(attr_name_);
  IODetails::DangerousString h5path(h5path_);
  parallelIO_write_simple_attr(attr_name.c_str(), (void*) &value,
          IODetails::PIO_DatatypeMap<scalar_type>::type, data_file_, h5path.c_str(), &IOgroup_);
}

template<typename scalar_type>
void
FileHDF5::WriteAttributeArray(scalar_type const * const values, std::size_t count,
        const std::string& attr_name, const std::string& h5path)
{
  Errors::Message message("FileHDF5::WriteAttributeArray not implemented for this type.");
  throw(message);
}


template<typename scalar_type>
scalar_type
FileHDF5::ReadAttribute(const std::string& attr_name_, const std::string& h5path_)
{
  IODetails::DangerousString attr_name(attr_name_);
  IODetails::DangerousString h5path(h5path_);

  scalar_type *loc_value;
  parallelIO_read_simple_attr(attr_name.c_str(), reinterpret_cast<void**>(&loc_value),
          IODetails::PIO_DatatypeMap<scalar_type>::type, data_file_, h5path.c_str(), &IOgroup_);
  scalar_type value = *loc_value;
  free(loc_value);
  return value;
}


template<typename scalar_type>
std::vector<scalar_type>
FileHDF5::ReadAttributeArray(std::size_t count,
                    const std::string& attr_name, const std::string& h5path)
{
  Errors::Message message("FileHDF5::ReadAttributeArray not implemented for this type.");
  throw(message);
}


template<typename scalar_type>
void
FileHDF5::WriteVector(const Vector_type_<scalar_type>& vec, const std::string& var_name)
{
  auto vec_data = vec.getData();

  IODetails::DangerousString full_h5path(var_name);

  // parallelIO only deals with global lengths that are ints
  // Really check this!
  static_assert(std::is_same<Tpetra::global_size_t, std::size_t>::value,
                "ASCEMIO does not support vectors with global length larger than a (signed) int");
  Tpetra::global_size_t globallength_raw = vec.getGlobalLength();

  if (globallength_raw > std::numeric_limits<int>::max()) {
    Errors::Message message("FileHDF5::ASCEMIO does not support vectors with global length larger than a (signed) int.");
    throw(message);
  }
  int globaldims[] = { static_cast<int>(globallength_raw), 1 };
  
  std::size_t locallength = vec.getLocalLength();
  int localdims[] = { static_cast<int>(locallength), 1 };
  
  parallelIO_write_dataset((void*) vec_data.get(),
                           IODetails::PIO_DatatypeMap<scalar_type>::type,
                           2, globaldims, localdims,
                           data_file_, full_h5path.c_str(),
                           &IOgroup_, NONUNIFORM_CONTIGUOUS_WRITE);
}

template<typename scalar_type>
void
FileHDF5::WriteMultiVectorBlock(const MultiVector_type_<scalar_type>& vec,
        const std::string& var_name)
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
  int globaldims[] = { static_cast<int>(vec.getNumVectors()), static_cast<int>(globallength_raw) };
  
  std::size_t locallength = vec.getLocalLength();
  int localdims[] = { globaldims[0], static_cast<int>(locallength) };

  IODetails::DangerousString full_h5path(var_name);
  auto vec_data = vec.get1dView();
  parallelIO_write_dataset((void*) vec_data.get(),
                           IODetails::PIO_DatatypeMap<scalar_type>::type,
                           2, globaldims, localdims,
                           data_file_, full_h5path.c_str(),
                           &IOgroup_, NONUNIFORM_CONTIGUOUS_WRITE);
}


template<typename scalar_type>
void
FileHDF5::WriteMultiVector(const MultiVector_type_<scalar_type>& vec,
                           const std::vector<std::string>& var_paths)
{
  if (var_paths.size() != vec.getNumVectors()) {
    Errors::Message message("FileHDF5::WriteMultiVector requested with invalid number of names.");
    throw(message);
  }

  // parallelIO only deals with global lengths that are ints
  // Really check this!
  static_assert(std::is_same<Tpetra::global_size_t, std::size_t>::value,
                "ASCEMIO does not support vectors with global length larger than a (signed) int");
  Tpetra::global_size_t globallength_raw = vec.getGlobalLength();

  if (globallength_raw > std::numeric_limits<int>::max()) {
    Errors::Message message("FileHDF5::ASCEMIO does not support vectors with global length larger than a (signed) int.");
    throw(message);
  }
  int globaldims[] = { static_cast<int>(globallength_raw), 1 };
  
  std::size_t locallength = vec.getLocalLength();
  int localdims[] = { static_cast<int>(locallength), 1 };

  for (std::size_t i=0; i!=vec.getNumVectors(); ++i) {
    IODetails::DangerousString full_h5path(var_paths[i]);

    auto vec_data = vec.getData(i);
    parallelIO_write_dataset((void*) vec_data.get(),
                             IODetails::PIO_DatatypeMap<scalar_type>::type,
                             2, globaldims, localdims,
                             data_file_, full_h5path.c_str(),
                             &IOgroup_, NONUNIFORM_CONTIGUOUS_WRITE);
  }
}


template<typename view_type>
void
FileHDF5::WriteView(const view_type vec, Tpetra::global_size_t global_length_raw,
        const std::string& var_path)
{
  Kokkos::View<double**, Kokkos::LayoutRight, AmanziDefaultHost> check_layout_right = vec;

  // parallelIO only deals with global lengths that are ints
  // Really check this!
  static_assert(std::is_same<Tpetra::global_size_t, std::size_t>::value,
                "ASCEMIO does not support vectors with global length larger than a (signed) int");
  if (global_length_raw > std::numeric_limits<int>::max()) {
    Errors::Message message("FileHDF5::ASCEMIO does not support vectors with global length larger than a (signed) int.");
    throw(message);
  }

  int globaldims[] = { static_cast<int>(global_length_raw), static_cast<int>(vec.extent(1)) };
  std::size_t locallength = vec.extent(0);
  int localdims[] = { static_cast<int>(locallength), globaldims[1] };

  IODetails::DangerousString full_h5path(var_path);
  parallelIO_write_dataset((void*) &vec(0,0),
                           PIO_DOUBLE,
                           2, globaldims, localdims,
                           data_file_, full_h5path.c_str(),
                           &IOgroup_, NONUNIFORM_CONTIGUOUS_WRITE);
}




template<typename scalar_type>
void
FileHDF5::ReadVector(Vector_type_<scalar_type>& vec,
        const std::string& var_name)
{
  IODetails::DangerousString full_h5path(var_name);

  int ndims;
  parallelIO_get_dataset_ndims(&ndims, data_file_, full_h5path.c_str(), &IOgroup_);
  if (ndims < 0) {
    Errors::Message message("FileHDF5::ReadVector data has negative dimension.");
    throw(message);
  }

  int globaldims[ndims], localdims[ndims];
  parallelIO_get_dataset_dims(globaldims, data_file_, full_h5path.c_str(), &IOgroup_);
  if (vec.getGlobalLength() != globaldims[0]) {
    Errors::Message message;
    message << "FileHDF5::ReadVector with incorrect length (got " << globaldims[0] << ", expected " << vec.getGlobalLength() << ").";
    throw(message);
  }

  localdims[0] = static_cast<int>(vec.getLocalLength());
  localdims[1] = globaldims[1];

  {
    auto vecv = vec.get1dViewNonConst();
    parallelIO_read_dataset((void *) vecv.get(),
                            IODetails::PIO_DatatypeMap<scalar_type>::type,
                            ndims, globaldims, localdims,
                            data_file_, full_h5path.c_str(),
                            &IOgroup_, NONUNIFORM_CONTIGUOUS_READ);
  }
}
    

template<typename scalar_type>
void
FileHDF5::ReadMultiVector(MultiVector_type_<scalar_type>& vec, const std::vector<std::string>& var_paths)
{
  {
    auto vecv = vec.get2dViewNonConst();

    if (var_paths.size() != vec.getNumVectors()) {
      Errors::Message message("FileHDF5::ReadMultiVector requested with invalid number of names.");
      throw(message);
    }
    
    for (std::size_t i=0; i!=vec.getNumVectors(); ++i) {
      IODetails::DangerousString full_h5path(var_paths[i]);

      int ndims;
      parallelIO_get_dataset_ndims(&ndims, data_file_, full_h5path.c_str(), &IOgroup_);
      if (ndims < 0) {
        Errors::Message message("FileHDF5::ReadVector data has negative dimension.");
        throw(message);
      }

      int globaldims[ndims], localdims[ndims];
      parallelIO_get_dataset_dims(globaldims, data_file_, full_h5path.c_str(), &IOgroup_);
      if (vec.getGlobalLength() != globaldims[0]) {
        Errors::Message message;
        message << "FileHDF5::ReadVector with incorrect length (got " << globaldims[0] << ", expected " << vec.getGlobalLength() << ").";
        throw(message);
      }

      localdims[0] = static_cast<int>(vec.getLocalLength());
      localdims[1] = globaldims[1];

      {
        parallelIO_read_dataset((void *) vecv[i].get(),
                IODetails::PIO_DatatypeMap<scalar_type>::type,
                ndims, globaldims, localdims,
                data_file_, full_h5path.c_str(),
                &IOgroup_, NONUNIFORM_CONTIGUOUS_READ);
      }
    }
  }
}



template<typename scalar_type>
void
FileHDF5::ReadMultiVectorBlock(MultiVector_type_<scalar_type>& vec, const std::string& var_name)
{
  {
    IODetails::DangerousString full_h5path(var_name);

    int ndims;
    parallelIO_get_dataset_ndims(&ndims, data_file_, full_h5path.c_str(), &IOgroup_);
    if (ndims < 0) {
      Errors::Message message("FileHDF5::ReadVector data has negative dimension.");
      throw(message);
    }

    int globaldims[ndims], localdims[ndims];
    parallelIO_get_dataset_dims(globaldims, data_file_, full_h5path.c_str(), &IOgroup_);

    // check dimensions
    if (vec.getGlobalLength() != globaldims[0]) {
      Errors::Message message;
      message << "FileHDF5::ReadVector with incorrect length (got " << globaldims[0] << ", expected " << vec.getGlobalLength() << ").";
      throw(message);
    }

    if (vec.getNumVectors() != globaldims[1]) {
      Errors::Message message;
      message << "FileHDF5::ReadVector with incorrect NumVectors (got " << globaldims[1] << ", expected " << vec.getNumVectors() << ").";
      throw(message);
    }

    localdims[0] = static_cast<int>(vec.getLocalLength());
    localdims[1] = globaldims[1];
    {
      auto vecv = vec.get1dViewNonConst();
      parallelIO_read_dataset((void *) vecv.get(),
              IODetails::PIO_DatatypeMap<scalar_type>::type,
              ndims, globaldims, localdims,
              data_file_, full_h5path.c_str(),
              &IOgroup_, NONUNIFORM_CONTIGUOUS_READ);
    }
  }
}


//
// STD::STRING
//
template<>
void
FileHDF5::WriteAttribute<std::string>(std::string value,
                             const std::string& attr_name, const std::string& h5path);
template<>
void
FileHDF5::WriteAttributeArray<std::string>(std::string const * const values, std::size_t count,
                          const std::string& attr_name, const std::string& h5path);
template<>
std::string
FileHDF5::ReadAttribute<std::string>(const std::string& attr_name, const std::string& h5path);
template<>
std::vector<std::string>
FileHDF5::ReadAttributeArray(std::size_t count,
                    const std::string& attr_name, const std::string& h5path);


  
} // namespace
