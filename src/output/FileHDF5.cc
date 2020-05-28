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

#include "dbc.hh"
#include "FileHDF5.hh"
#include "AmanziVector.hh"

namespace Amanzi {

FileHDF5::FileHDF5(const Comm_ptr_type& comm, const std::string& filename,
                   file_mode_t mode)
  : comm_(comm), filename_(filename), data_file_(-1)
{
  auto mpi_comm = Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm_);
  AMANZI_ASSERT(mpi_comm.get());

  IOconfig_.numIOgroups = 1;
  IOconfig_.commIncoming = *mpi_comm->getRawMpiComm();
  parallelIO_IOgroup_init(&IOconfig_, &IOgroup_);
  OpenFile(mode);
}

FileHDF5::~FileHDF5()
{
  CloseFile();
  parallelIO_IOgroup_cleanup(&IOgroup_);
}


void
FileHDF5::OpenFile(file_mode_t mode)
{
  if (data_file_ < 0) {
    data_file_ = parallelIO_open_file(filename_.c_str(), &IOgroup_, mode);
    if (data_file_ < 0) {
      Errors::Message message;
      message << "FileHDF5: Unable to open file \"" << filename_
              << "\" with HDF5 error code " << (int) data_file_;
      throw(message);
    }
  }
}

void
FileHDF5::CloseFile()
{
  if (open_groups_.size() > 0) {
    AMANZI_ASSERT(data_file_ >= 0);
  }
  for (auto gid : open_groups_) {
    int ierr = parallelIO_close_dataset_group(data_file_, &IOgroup_);
    AMANZI_ASSERT(!ierr);
  }
  open_groups_.clear();

  if (data_file_ >= 0) {
    int ierr = parallelIO_close_file(data_file_, &IOgroup_);
    AMANZI_ASSERT(!ierr);
    data_file_ = -1;
  }
}

void
FileHDF5::CreateGroup(const std::string& h5path_)
{
  IODetails::DangerousString h5path(h5path_);
  if (open_groups_.size() > 0) {
    // NOTE: this is a fairly straightforward fix of ASCEMIO, whose data
    // structures currently only store one groupid.  It would be
    // straightforward to make that an array (much like file ids are an array),
    // but it is not done. --etc
    Errors::Message message;
    message << "ASCEMIO currently limited to 1 open group per file at a time: \"" << filename_ << "\"";
    throw(message);
  }
  open_groups_.push_back(
    parallelIO_create_dataset_group(h5path.c_str(), data_file_, &IOgroup_));
  all_groups_.push_back(h5path_);
}

void
FileHDF5::CloseGroup()
{
  if (open_groups_.size() == 1) {
    parallelIO_close_dataset_group(data_file_, &IOgroup_);
    open_groups_.clear();
  }
}

bool
FileHDF5::HasGroup(const std::string& h5path) const
{
  return std::find(all_groups_.begin(), all_groups_.end(), h5path) != all_groups_.end();
}
  

void
FileHDF5::WriteVector(const std::string& var_name, const Map_type& vec)
{
  IODetails::DangerousString full_h5path(var_name);

  // parallelIO only deals with global lengths that are ints
  // Really check this!
  static_assert(std::is_same<Tpetra::global_size_t, std::size_t>::value,
                "ASCEMIO does not support vectors with global length larger "
                "than a (signed) int");
  auto globallength_raw = vec.getGlobalNumElements();

  if (globallength_raw > std::numeric_limits<int>::max()) {
    Errors::Message message("FileHDF5::ASCEMIO does not support vectors with "
                            "global length larger than a (signed) int.");
    throw(message);
  }
  int globaldims[] = { static_cast<int>(globallength_raw), 1 };

  std::size_t locallength = vec.getNodeNumElements();
  int localdims[] = { static_cast<int>(locallength), 1 };

  Teuchos::Array<int> map_data(vec.getNodeNumElements());
  {
    auto vec_data = vec.getMyGlobalIndices();
    for (int i = 0; i != map_data.size(); ++i) { map_data[i] = vec_data[i]; }
  }
  parallelIO_write_dataset((void*)map_data.data(),
                           PIO_INTEGER,
                           2,
                           globaldims,
                           localdims,
                           data_file_,
                           full_h5path.c_str(),
                           &IOgroup_,
                           NONUNIFORM_CONTIGUOUS_WRITE);
}


//
// STD::STRING
//
template <>
void
FileHDF5::WriteAttribute<std::string>(const std::string& attr_name_,
                                      const std::string& h5path_,
                                      std::string value)
{
  if (!HasGroup(h5path_)) {
    CreateGroup(h5path_);
    CloseGroup();
  }
  
  IODetails::DangerousString attr_name(attr_name_);
  IODetails::DangerousString h5path(h5path_);
  parallelIO_write_simple_attr(attr_name.c_str(),
                               (void*)value.c_str(),
                               PIO_STRING,
                               data_file_,
                               h5path.c_str(),
                               &IOgroup_);
}


template <>
std::string
FileHDF5::ReadAttribute<std::string>(const std::string& attr_name_,
                                     const std::string& h5path_)
{
  IODetails::DangerousString attr_name(attr_name_);
  IODetails::DangerousString h5path(h5path_);

  char* loc_value;
  parallelIO_read_simple_attr(attr_name.c_str(),
                              reinterpret_cast<void**>(&loc_value),
                              PIO_STRING,
                              data_file_,
                              h5path.c_str(),
                              &IOgroup_);
  std::string value(loc_value);
  free(loc_value);
  return value;
}

} // namespace Amanzi
