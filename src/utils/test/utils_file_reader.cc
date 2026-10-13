/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Tests FileReader (FileReaderHDF5/FileReaderNetCDF, via createFileReader):
//! round-trips a small 1D array and a small paged 3D chunk through each
//! backend, confirming both agree.

#include "UnitTest++.h"

#include <string>

#include "AmanziComm.hh"
#include "FileHDF5.hh"
#include "FileNetCDF.hh"
#include "FileReader.hh"
#include "IOFileFormat.hh"

using namespace Amanzi;
using namespace Amanzi::Utils;

void
run_read1D_test(IOFileFormat_kind file_format, const std::string& filename)
{
  // Write a small 1D array directly through the backend, then read it
  // back through the FileReader interface.
  Kokkos::View<double*, Kokkos::HostSpace> written("written", 5);
  for (int i = 0; i != 5; ++i) written(i) = 10. * i + 1.;

  if (file_format == IOFileFormat_kind::HDF5) {
    FileHDF5 out(getDefaultComm(), filename, FILE_CREATE);
    out.writeLocalView("x", "/", written);
  } else {
    FileNetCDF out(getDefaultComm(), filename, FILE_CREATE);
    out.writeLocalView("x", written);
  }

  auto reader = createFileReader(file_format, getDefaultComm(), filename);
  Kokkos::View<double*, Kokkos::HostSpace> read_back;
  reader->read1D("x", read_back);

  CHECK_EQUAL(5, (int)read_back.extent(0));
  for (int i = 0; i != 5; ++i) CHECK_CLOSE(written(i), read_back(i), 1.e-12);
}


TEST(file_reader_read1d_hdf5)
{
  run_read1D_test(IOFileFormat_kind::HDF5, "test_file_reader_1d.h5");
}


TEST(file_reader_read1d_netcdf)
{
  run_read1D_test(IOFileFormat_kind::NETCDF, "test_file_reader_1d.nc");
}


void
run_paged_chunk_test(IOFileFormat_kind file_format, const std::string& filename)
{
  // Write a small [NTIMES, NX, NY] array directly through the backend,
  // then read a sub-chunk of its time axis back through FileReader.
  const int ntimes = 4, nx = 2, ny = 3;
  Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> written(
    "written", ntimes, nx, ny);
  for (int t = 0; t != ntimes; ++t)
    for (int i = 0; i != nx; ++i)
      for (int j = 0; j != ny; ++j) written(t, i, j) = 100. * t + 10. * i + j;

  if (file_format == IOFileFormat_kind::HDF5) {
    FileHDF5 out(getDefaultComm(), filename, FILE_CREATE);
    out.writeLocalView("v", "/", written, FileHDF5::TimeRange{ 0, (std::size_t)ntimes });
  } else {
    FileNetCDF out(getDefaultComm(), filename, FILE_CREATE);
    out.writeLocalView("v", written, FileNetCDF::TimeRange{ 0, (std::size_t)ntimes });
  }

  auto reader = createFileReader(file_format, getDefaultComm(), filename);
  Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> chunk("chunk", 2, nx, ny);
  reader->readPagedChunk("v", FileReader::TimeRange{ 1, 3 }, chunk);

  for (int t = 0; t != 2; ++t)
    for (int i = 0; i != nx; ++i)
      for (int j = 0; j != ny; ++j) CHECK_CLOSE(written(t + 1, i, j), chunk(t, i, j), 1.e-12);
}


TEST(file_reader_paged_chunk_hdf5)
{
  run_paged_chunk_test(IOFileFormat_kind::HDF5, "test_file_reader_chunk.h5");
}


TEST(file_reader_paged_chunk_netcdf)
{
  run_paged_chunk_test(IOFileFormat_kind::NETCDF, "test_file_reader_chunk.nc");
}
