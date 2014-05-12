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
  if (op_list.isSublist("diffusion operator")) {
    Teuchos::ParameterList dlist = op_list.sublist("diffusion operator");

    std::vector<std::string> names;
    names = dlist.get<Teuchos::Array<std::string> > ("schema").toVector();

    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh);
    cvs->SetGhosted(true);

    if (names[0] == "cell") {
       cvs->SetComponent("cell", AmanziMesh::CELL, 1);
    } else if (names[0] == "node") {
       cvs->SetComponent("node", AmanziMesh::NODE, 1);
    }

    for (int i = 1; i < names.size(); i++) {
      cvs->SetOwned(false);
      if (names[i] == "cell") {
         cvs->AddComponent("cell", AmanziMesh::CELL, 1);
      } else if (names[i] == "node") {
         cvs->AddComponent("node", AmanziMesh::NODE, 1);
      } else if (names[i] == "face") {
         cvs->AddComponent("face", AmanziMesh::FACE, 1);
      }
    }

    Teuchos::RCP<OperatorDiffusion> op = Teuchos::rcp(new OperatorDiffusion(cvs, dlist));
    op->Init();
    return op;

  } else {
    Errors::Message msg("OperatorDiffusionFactory: \"diffusion operator\" does not exist.");
    Exceptions::amanzi_throw(msg);
  }
}

}  // namespace Operators
}  // namespace Amanzi


