//! Factory for Input objects.
/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

/*

*/

#include <string>

#include "errors.hh"
#include "Input.hh"
#include "InputOutputHDF5.hh"

#include "InputFactory.hh"

namespace Amanzi {

Teuchos::RCP<Input>
CreateInput(Teuchos::ParameterList& plist,
        const Comm_ptr_type& comm)
{
  std::string input_type = plist.get<std::string>("input type");
  if (input_type == "HDF5") {
    return Teuchos::rcp(new InputOutputHDF5(plist, comm));
  } else {
    Errors::Message msg;
    msg << "InputFactory: Unknown input type \"" << input_type << "\"";
    throw(msg);
  }
}


} // namespace
