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

#ifndef AMANZI_INPUT_FACTORY_HH_
#define AMANZI_INPUT_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AmanziTypes.hh"

namespace Amanzi {

class Input;

namespace AmanziMesh {
class Mesh;
}

namespace InputFactory {

std::unique_ptr<Input>
createForCheckpoint(Teuchos::ParameterList& plist, const Comm_ptr_type& comm);

// Dispatches on plist's "file format" ("HDF5", the default, or "NetCDF"),
// constructing an InputHDF5 or InputVerticallyStructuredNetCDF as
// appropriate.  mesh is only used by the NetCDF case, to derive columnar
// mesh structure.
std::unique_ptr<Input>
createForFile(Teuchos::ParameterList& plist,
              const Comm_ptr_type& comm,
              const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

} // namespace InputFactory
} // namespace Amanzi

#endif
