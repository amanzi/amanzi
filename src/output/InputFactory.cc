/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Factory for Input objects.
/*

*/

#include <string>

#include "errors.hh"
#include "IOFileFormat.hh"
#include "Input.hh"
#include "InputOutputHDF5.hh"
#include "InputOutputVerticallyStructuredNetCDF.hh"

#include "InputFactory.hh"

namespace Amanzi {
namespace InputFactory {

std::unique_ptr<Input>
createForCheckpoint(Teuchos::ParameterList& plist, const Comm_ptr_type& comm)
{
  std::string filename = plist.get<std::string>("file name");
  IOFileFormat_kind file_format = createIOFileFormat(plist.get<std::string>("file format", "HDF5"));

  switch (file_format) {
  case IOFileFormat_kind::HDF5:
    return std::make_unique<InputHDF5>(comm, filename);
  default: {
    Errors::Message msg;
    msg << "InputFactory::createForCheckpoint: Unknown input type \"" << to_string(file_format)
        << "\"";
    throw(msg);
  }
  }
}


std::unique_ptr<Input>
createForFile(Teuchos::ParameterList& plist,
              const Comm_ptr_type& comm,
              const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
{
  std::string filename = plist.get<std::string>("file name");
  IOFileFormat_kind file_format = createIOFileFormat(plist.get<std::string>("file format", "HDF5"));

  switch (file_format) {
  case IOFileFormat_kind::HDF5:
    return std::make_unique<InputHDF5>(comm, filename);
  case IOFileFormat_kind::NETCDF:
    return std::make_unique<InputVerticallyStructuredNetCDF>(comm, filename, mesh);
  default: {
    Errors::Message msg;
    msg << "InputFactory::createForFile: Unknown input type \"" << to_string(file_format) << "\"";
    throw(msg);
  }
  }
}

} // namespace InputFactory
} // namespace Amanzi
