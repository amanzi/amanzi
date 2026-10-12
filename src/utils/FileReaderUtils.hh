/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! readLocalArray: one-shot read of a 1D/2D array from an HDF5 or NetCDF file.
/*
  For callers that construct a reader, read once or twice, and let it go
  out of scope -- no held reader, no paging.
*/
#pragma once

#include <string>

#include "AmanziTypes.hh"
#include "IOFileFormat.hh"
#include "FileHDF5.hh"
#include "FileNetCDF.hh"

namespace Amanzi {

template <typename view_type>
void
readLocalArray(IOFileFormat_kind file_format,
               const Comm_ptr_type& comm,
               const std::string& filename,
               const std::string& name,
               view_type& view)
{
  switch (file_format) {
  case IOFileFormat_kind::HDF5: {
    FileHDF5 file(comm, filename, FILE_READONLY);
    file.readLocalView(name, "/", view);
    break;
  }
  case IOFileFormat_kind::NETCDF: {
    FileNetCDF file(comm, filename, FILE_READONLY);
    file.readLocalView(name, view);
    break;
  }
  }
}

} // namespace Amanzi
