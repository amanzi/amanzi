/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  MPC PK

  Weak coupling of mechanics and coupled flow PKs.
*/

#include <string>

#include "PK_Utils.hh"
#include "Transport_PK.hh"

#include "MechanicsFlowMatrixFracture_PK.hh"

namespace Amanzi {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* Constructor
****************************************************************** */
MechanicsFlowMatrixFracture_PK::MechanicsFlowMatrixFracture_PK(
  Teuchos::ParameterList& pk_tree,
  const Teuchos::RCP<Teuchos::ParameterList>& global_list,
  const Teuchos::RCP<State>& S,
  const Teuchos::RCP<TreeVector>& soln)
  : PK_MPCWeak(pk_tree, global_list, S, soln)
{
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = my_list_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("MechanicsCoupledFlow", vlist));
}


/* ******************************************************************
* Setup of PK
****************************************************************** */
void
MechanicsFlowMatrixFracture_PK::Setup()
{
  PK_MPCWeak::Setup();

  std::string passwd("");
  hydrostatic_stress_key_ = "hydrostatic_stress";

  S_->Require<CV_t, CVS_t>(hydrostatic_stress_key_, Tags::DEFAULT, passwd)
    .SetMesh(S_->GetMesh())
    ->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
}

} // namespace Amanzi
