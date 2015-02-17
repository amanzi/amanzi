/*
  This is the operators component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base factory for diffusion operators.
*/

#include "errors.hh"

#include "BCs.hh"
#include "OperatorDefs.hh"
#include "OperatorDiffusionFactory.hh"
#include "OperatorDiffusion.hh"
//#include "OperatorDiffusionTPFA.hh"
#include "OperatorDiffusionWithGravity.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
 * Initialization of the diffusion operators.
 ****************************************************************** */
Teuchos::RCP<OperatorDiffusion> OperatorDiffusionFactory::Create(
    Teuchos::RCP<const AmanziMesh::Mesh> mesh, 
    Teuchos::RCP<BCs> bc, 
    Teuchos::ParameterList& oplist,
    const AmanziGeometry::Point& g,
    int upwind_method)
{
  std::vector<std::string> names;
  names = oplist.get<Teuchos::Array<std::string> > ("schema").toVector();
  int nnames = names.size();

  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh);
  cvs->SetGhosted(true);

  std::vector<AmanziMesh::Entity_kind> locations(nnames);
  std::vector<int> num_dofs(nnames, 1);
 
  for (int i = 0; i < nnames; i++) {
    if (names[i] == "cell") {
      locations[i] = AmanziMesh::CELL;
    } else if (names[i] == "node") {
      locations[i] = AmanziMesh::NODE;
    } else if (names[i] == "face") {
      locations[i] = AmanziMesh::FACE;
    }
  }

  cvs->SetComponents(names, locations, num_dofs);
  cvs->SetOwned(false);

  // Let us try to identify a FV scheme.
  std::string name = oplist.get<std::string>("discretization primary");
  if (name == "fv: default") {

    Errors::Message msg("OperatorDiffusionFactory: TPFA not yet implemented");
    Exceptions::amanzi_throw(msg);
  }

  // Let us see if we have gravity.
  bool flag = oplist.get<bool>("gravity", false);
  if (!flag) {
    return Teuchos::rcp(new OperatorDiffusion(oplist, mesh));
  } else {
    Teuchos::RCP<OperatorDiffusionWithGravity> op = Teuchos::rcp(new OperatorDiffusionWithGravity(oplist, mesh));
    op->SetGravity(g);
    return op;
  }
}

}  // namespace Operators
}  // namespace Amanzi


