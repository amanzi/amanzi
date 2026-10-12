/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! FileReader: holds either a FileHDF5 or a FileNetCDF, chosen at runtime.
/*
  For a caller that keeps a reader alive across repeated calls, as
  opposed to readLocalArray (FileReaderUtils.hh), which is for one-shot
  construct-read-destroy callers.
*/
#pragma once

#include <memory>
#include <string>

#include "Kokkos_Core.hpp"

#include "AmanziTypes.hh"
#include "IOFileFormat.hh"

namespace Amanzi {

class FileReader {
 public:
  // [begin, end) range of time indices, for readPagedChunk.
  using TimeRange = Kokkos::pair<std::size_t, std::size_t>;

  virtual ~FileReader() = default;

  // Unpaged, read once (e.g. the x/y/time grid).  view may be
  // default-constructed (size 0); it is then resized to fit the on-disk
  // shape.
  virtual void
  read1D(const std::string& name, Kokkos::View<double*, Kokkos::HostSpace>& view) const = 0;

  // Paged: on-disk array has one leading time dimension beyond view's own
  // (2D) shape; view is rank 3, dimension 0 selecting time_chunk.
  virtual void
  readPagedChunk(const std::string& name,
                 TimeRange time_chunk,
                 Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace>& view) const = 0;
};

std::unique_ptr<FileReader> createFileReader(IOFileFormat_kind file_format,
                                             const Comm_ptr_type& comm,
                                             const std::string& filename);

} // namespace Amanzi
