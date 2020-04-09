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

#ifndef AMANZI_OUTPUT_FACTORY_HH_
#define AMANZI_OUTPUT_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AmanziTypes.hh"

namespace Amanzi {

class Output;
namespace AmanziMesh {
class Mesh;
}

namespace OutputFactory {


std::unique_ptr<Output>
CreateForVis(Teuchos::ParameterList& plist,
             const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

std::unique_ptr<Output>
CreateForCheckpoint(Teuchos::ParameterList& plist, const Comm_ptr_type& comm);

} // namespace OutputFactory
} // namespace Amanzi

#endif
