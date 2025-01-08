/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "errors.hh"
#include "Key.hh"
#include "Reader.hh"
#include "ReaderHDF5.hh"
#include "ReaderNetCDF.hh"

namespace Amanzi {

std::unique_ptr<Reader>
createReader(const std::string& filename)
{
  if (Keys::ends_with(filename, ".h5")) {
    return std::make_unique<ReaderHDF5>(filename);
  } else if (Keys::ends_with(filename, ".nc")) {
    return std::make_unique<ReaderNetCDF>(filename);
  } else {
    Errors::Message msg;
    msg << "Unable to create a reader for \"" << filename
        << "\" -- valid extensions are .nc and .h5";
    Exceptions::amanzi_throw(msg);
  }
  return nullptr;
}

} // namespace Amanzi
