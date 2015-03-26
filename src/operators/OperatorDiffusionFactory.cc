/*
  This is the operators component of the Amanzi code.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  Base factory for diffusion operators.
*/

#include "errors.hh"

#include "BCs.hh"
#include "OperatorDefs.hh"
#include "OperatorDiffusionFactory.hh"
#include "OperatorDiffusionMFD.hh"
#include "OperatorDiffusionFV.hh"
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
  // Let us try to identify a FV scheme.
  std::string name = oplist.get<std::string>("discretization primary");
  if (name == "fv: default") {
    Teuchos::RCP<OperatorDiffusionFV> op = Teuchos::rcp(new OperatorDiffusionFV(oplist, mesh));
    if (oplist.get<bool>("gravity", false)) op->SetGravity(g);
    op->SetBCs(bc);
    return op;
  }

  // Let us see if we have gravity.
  bool flag = oplist.get<bool>("gravity", false);
  if (!flag) {
    Teuchos::RCP<OperatorDiffusionMFD> op = Teuchos::rcp(new OperatorDiffusionMFD(oplist, mesh));
    op->SetBCs(bc);
    return op;
  } else {
    Teuchos::RCP<OperatorDiffusionWithGravity> op = Teuchos::rcp(new OperatorDiffusionWithGravity(oplist, mesh));
    op->SetGravity(g);
    op->SetBCs(bc);
    return op;
  }
}

}  // namespace Operators
}  // namespace Amanzi


