/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Factory for Input objects.

/*

*/

#include <string>

#include "errors.hh"
#include "UniqueHelpers.hh"
#include "Input.hh"
#include "InputOutputHDF5.hh"

#include "InputFactory.hh"

namespace Amanzi {
namespace InputFactory {

std::unique_ptr<Input>
CreateForCheckpoint(Teuchos::ParameterList& plist, const Comm_ptr_type& comm)
{
  std::string filename = plist.get<std::string>("file name");
  std::string input_type = plist.get<std::string>("file format", "HDF5");

  if (input_type == "HDF5") {
    return std::make_unique<InputHDF5>(comm, filename);
  } else {
    Errors::Message msg;
    msg << "InputFactory: Unknown input type \"" << input_type << "\"";
    throw(msg);
  }
}

} // namespace InputFactory
} // namespace Amanzi
