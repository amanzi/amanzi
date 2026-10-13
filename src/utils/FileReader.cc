/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "errors.hh"
#include "exceptions.hh"
#include "FileHDF5.hh"
#include "FileNetCDF.hh"
#include "FileReader.hh"

namespace Amanzi {
namespace Utils {

namespace {

class FileReaderHDF5 : public FileReader {
 public:
  FileReaderHDF5(const Comm_ptr_type& comm, const std::string& filename)
    : file_(comm, filename, FILE_READONLY)
  {}

  void
  read1D(const std::string& name, Kokkos::View<double*, Kokkos::HostSpace>& view) const override
  {
    file_.readLocalView(name, "/", view);
  }

  void readPagedChunk(
    const std::string& name,
    TimeRange time_chunk,
    Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace>& view) const override
  {
    file_.readLocalView(
      name, "/", view, FileHDF5::TimeRange{ time_chunk.first, time_chunk.second });
  }

 private:
  mutable FileHDF5 file_;
};

class FileReaderNetCDF : public FileReader {
 public:
  FileReaderNetCDF(const Comm_ptr_type& comm, const std::string& filename)
    : file_(comm, filename, FILE_READONLY)
  {}

  void
  read1D(const std::string& name, Kokkos::View<double*, Kokkos::HostSpace>& view) const override
  {
    file_.readLocalView(name, view);
  }

  void readPagedChunk(
    const std::string& name,
    TimeRange time_chunk,
    Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace>& view) const override
  {
    file_.readLocalView(name, view, FileNetCDF::TimeRange{ time_chunk.first, time_chunk.second });
  }

 private:
  mutable FileNetCDF file_;
};

} // namespace

std::unique_ptr<FileReader>
createFileReader(IOFileFormat_kind file_format,
                 const Comm_ptr_type& comm,
                 const std::string& filename)
{
  switch (file_format) {
  case IOFileFormat_kind::HDF5:
    return std::make_unique<FileReaderHDF5>(comm, filename);
  case IOFileFormat_kind::NETCDF:
    return std::make_unique<FileReaderNetCDF>(comm, filename);
  default: {
    Errors::Message msg;
    msg << "createFileReader: invalid \"file format\" \"" << to_string(file_format)
        << "\" for reading a file";
    Exceptions::amanzi_throw(msg);
    return nullptr;
  }
  }
}

} // namespace Utils
} // namespace Amanzi
