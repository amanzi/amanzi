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

#ifndef AMANZI_OUTPUT_FACTORY_HH_
#define AMANZI_OUTPUT_FACTORY_HH_

#include "Teuchos_RCP.hpp"

namespace Amanzi {

class Output;
namespace AmanziMesh {
class Mesh;
}

Teuchos::RCP<Output> CreateOutput(Teuchos::ParameterList& plist,
        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);


} // namespace

#endif
