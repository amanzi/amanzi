/*
  This is the operators component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base factory for diffusion operators.
  Usage:
*/

#include "errors.hh"

#include "OperatorDefs.hh"
#include "OperatorDiffusionFactory.hh"
#include "OperatorDiffusion.hh"
#include "OperatorDiffusionSurface.hh"


namespace Amanzi {
namespace Operators {

/* ******************************************************************
 * Initialization of the diffusion operators.
 ****************************************************************** */
Teuchos::RCP<OperatorDiffusion> OperatorDiffusionFactory::Create(
    Teuchos::RCP<const AmanziMesh::Mesh> mesh, const Teuchos::ParameterList& op_list)
{
  if (op_list.isParameter("operator type")) {
    std::string type = op_list.get<std::string>("operator type");

    if (type == "mfd nodal domain") {
      Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
      cvs->SetMesh(mesh);
      cvs->SetGhosted(true);
      cvs->SetComponent("node", AmanziMesh::NODE, 1);

      Teuchos::RCP<OperatorDiffusion> op = Teuchos::rcp(new OperatorDiffusion(cvs, 0));
      op->Init();
      return op;
    }
  } else {
    Errors::Message msg("OperatorDiffusionFactory: wrong value of parameter `\"operator type`\"");
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace Operators
}  // namespace Amanzi


