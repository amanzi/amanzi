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
  Utils::IOFileFormat_kind file_format = Utils::createIOFileFormat(plist, filename);

  switch (file_format) {
  case Utils::IOFileFormat_kind::HDF5:
    return std::make_unique<InputHDF5>(comm, filename);
  default: {
    Errors::Message msg;
    msg << "InputFactory::createForCheckpoint: invalid \"file format\" \"" << to_string(file_format)
        << "\" for reading a checkpoint file";
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
  Utils::IOFileFormat_kind file_format = Utils::createIOFileFormat(plist, filename);

  switch (file_format) {
  case Utils::IOFileFormat_kind::HDF5:
    return std::make_unique<InputHDF5>(comm, filename);
  case Utils::IOFileFormat_kind::NETCDF:
    return std::make_unique<InputVerticallyStructuredNetCDF>(comm, filename, mesh);
  default: {
    Errors::Message msg;
    msg << "InputFactory::createForFile: invalid \"file format\" \"" << to_string(file_format)
        << "\" for reading a file";
    throw(msg);
  }
  }
}

} // namespace InputFactory
} // namespace Amanzi
