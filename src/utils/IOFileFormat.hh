/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! IOFileFormat_kind: which low-level file backend (HDF5 or NetCDF) to use.
#pragma once

#include <string>

#include "errors.hh"
#include "exceptions.hh"

namespace Amanzi {

enum class IOFileFormat_kind { HDF5, NETCDF };

inline std::string
to_string(const IOFileFormat_kind format)
{
  switch (format) {
  case (IOFileFormat_kind::HDF5):
    return "HDF5";
  case (IOFileFormat_kind::NETCDF):
    return "NetCDF";
  default:
    return "unknown";
  }
}

inline IOFileFormat_kind
createIOFileFormat(const std::string& format_string)
{
  if (format_string == "HDF5") {
    return IOFileFormat_kind::HDF5;
  } else if (format_string == "NetCDF") {
    return IOFileFormat_kind::NETCDF;
  } else {
    Errors::Message msg;
    msg << "Unknown \"file format\" string: \"" << format_string
        << "\", valid are \"HDF5\" and \"NetCDF\".";
    Exceptions::amanzi_throw(msg);
    return IOFileFormat_kind::HDF5;
  }
}

} // namespace Amanzi
