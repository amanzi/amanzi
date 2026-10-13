/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! IOFileFormat_kind: which file format/backend to use for a given file.
#pragma once

#include <algorithm>
#include <string>

#include "Teuchos_ParameterList.hpp"

#include "errors.hh"
#include "exceptions.hh"

namespace Amanzi {
namespace Utils {

enum class IOFileFormat_kind { HDF5, NETCDF, XDMF, SILO };

inline std::string
to_string(const IOFileFormat_kind format)
{
  switch (format) {
  case (IOFileFormat_kind::HDF5):
    return "HDF5";
  case (IOFileFormat_kind::NETCDF):
    return "NetCDF";
  case (IOFileFormat_kind::XDMF):
    return "XDMF";
  case (IOFileFormat_kind::SILO):
    return "SILO";
  default:
    return "unknown";
  }
}

inline IOFileFormat_kind
createIOFileFormatFromString(const std::string& format_string)
{
  std::string upper(format_string);
  std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
  if (upper == "HDF5") {
    return IOFileFormat_kind::HDF5;
  } else if (upper == "NETCDF") {
    return IOFileFormat_kind::NETCDF;
  } else if (upper == "XDMF") {
    return IOFileFormat_kind::XDMF;
  } else if (upper == "SILO") {
    return IOFileFormat_kind::SILO;
  } else {
    Errors::Message msg;
    msg << "Unknown \"file format\" string: \"" << format_string
        << "\", valid are \"HDF5\", \"NetCDF\", \"XDMF\", and \"SILO\".";
    Exceptions::amanzi_throw(msg);
    return IOFileFormat_kind::HDF5;
  }
}

// Infers the file format from the filename's extension (.h5 --> HDF5, .nc
// --> NetCDF, .xmf --> XDMF, .silo --> SILO).  If the extension is not
// recognized, "file format" must be provided in plist.
inline IOFileFormat_kind
createIOFileFormat(const Teuchos::ParameterList& plist, const std::string& filename)
{
  auto ext_pos = filename.rfind('.');
  std::string ext = ext_pos == std::string::npos ? "" : filename.substr(ext_pos);
  if (ext == ".h5") {
    return IOFileFormat_kind::HDF5;
  } else if (ext == ".nc") {
    return IOFileFormat_kind::NETCDF;
  } else if (ext == ".xmf") {
    return IOFileFormat_kind::XDMF;
  } else if (ext == ".silo") {
    return IOFileFormat_kind::SILO;
  } else if (plist.isParameter("file format")) {
    return createIOFileFormatFromString(plist.get<std::string>("file format"));
  } else {
    Errors::Message msg;
    msg << "Cannot infer \"file format\" from filename \"" << filename
        << "\" -- must be provided explicitly as \"file format\".";
    Exceptions::amanzi_throw(msg);
    return IOFileFormat_kind::HDF5;
  }
}

// For output, where no on-disk filename exists yet to infer from -- reads
// "file format" from plist, falling back to default_format if absent.
inline IOFileFormat_kind
createIOFileFormat(const Teuchos::ParameterList& plist, IOFileFormat_kind default_format)
{
  if (plist.isParameter("file format")) {
    return createIOFileFormatFromString(plist.get<std::string>("file format"));
  }
  return default_format;
}

} // namespace Utils
} // namespace Amanzi
