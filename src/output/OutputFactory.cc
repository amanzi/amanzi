//! Factory for Output objects.
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
#include "Output.hh"
#include "OutputXDMF.hh"
#include "InputOutputHDF5.hh"

#include "OutputFactory.hh"

namespace Amanzi {

Teuchos::RCP<Output>
CreateOutput(Teuchos::ParameterList& plist,
        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
{
  std::string output_type = plist.get<std::string>("output type");
  if (output_type == "XDMF") {
    return Teuchos::rcp(new OutputXDMF(plist, mesh));
  } else if (output_type == "HDF5") {
    return Teuchos::rcp(new InputOutputHDF5(plist, mesh->get_comm()));
  } else {
    Errors::Message msg;
    msg << "OutputFactory: Unknown output type \"" << output_type << "\"";
    throw(msg);
  }
}


} // namespace
