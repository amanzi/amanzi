/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Factory for Output objects.

/*

*/

#include <string>

#include "errors.hh"
#include "Output.hh"
#include "OutputXDMF.hh"
#include "InputOutputHDF5.hh"

#include "OutputFactory.hh"


namespace Amanzi {
namespace OutputFactory {

std::unique_ptr<Output>
CreateForVis(Teuchos::ParameterList& plist,
             const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
{
  std::string output_type = plist.get<std::string>("file format", "XDMF");
  if (output_type == "XDMF") {
    return std::make_unique<OutputXDMF>(plist, mesh);
  } else {
    Errors::Message msg;
    msg << "OutputFactory: Unknown output type \"" << output_type << "\"";
    throw(msg);
  }
}

std::unique_ptr<Output>
CreateForCheckpoint(Teuchos::ParameterList& plist, const Comm_ptr_type& comm)
{
  std::string output_type = plist.get<std::string>("file format", "HDF5");
  if (output_type == "HDF5") {
    return std::make_unique<OutputHDF5>(plist, comm);
  } else {
    Errors::Message msg;
    msg << "OutputFactory: Unknown output type \"" << output_type << "\"";
    throw(msg);
  }
}


} // namespace OutputFactory
} // namespace Amanzi
