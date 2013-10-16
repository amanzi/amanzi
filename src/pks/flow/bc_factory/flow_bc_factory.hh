#ifndef AMANZI_FLOW_BC_FACTORY_HH_
#define AMANZI_FLOW_BC_FACTORY_HH_

/* -------------------------------------------------------------------------
ATS

Author: ...
    Ethan Coon (ATS version) (ecoon@lanl.gov)

*/

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "bc_factory.hh"

namespace Amanzi {
namespace Flow {

class FlowBCFactory : public Amanzi::BCFactory {

public:
  FlowBCFactory(const Teuchos::RCP<const AmanziMesh::Mesh> &mesh,
                const Teuchos::ParameterList& plist)
      : Amanzi::BCFactory(mesh,plist) {}

  Teuchos::RCP<Functions::BoundaryFunction> CreatePressure() const {
    return CreateWithFunction("pressure", "boundary pressure");
  }

  Teuchos::RCP<Functions::BoundaryFunction> CreateHead() const {
    return CreateWithFunction("head", "boundary head");
  }

  Teuchos::RCP<Functions::BoundaryFunction> CreateMassFlux() const {
    return CreateWithFunction("mass flux", "outward mass flux");
  }

  Teuchos::RCP<Functions::BoundaryFunction> CreateZeroGradient() const {
    return CreateWithoutFunction("zero gradient");
  }

  Teuchos::RCP<Functions::BoundaryFunction> CreateSeepageFace() const {
    return CreateWithFunction("seepage face", "outward mass flux");
  }

};

}  // namespace
}  // namespace

#endif // AMANZI_FLOW_BC_FACTORY_HH_
