/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_FLOW_DOMAIN_FUNCTION_HH_
#define AMANZI_FLOW_DOMAIN_FUNCTION_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "Mesh.hh"
#include "PK_DomainFunction.hh"

namespace Amanzi {
namespace Flow {

class FlowDomainFunction : public PK_DomainFunction {
 public:
  FlowDomainFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      PK_DomainFunction(mesh) {};
  ~FlowDomainFunction() {};
};

}  // namespace Flow
}  // namespace Amanzi

#endif

