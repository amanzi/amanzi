/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Mechanics PK

*/

#include <vector>

#include "MechanicsFracturedMatrix_PK.hh"

namespace Amanzi {
namespace Mechanics {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* New constructor: extracts lists and requires fields.
****************************************************************** */
MechanicsFracturedMatrix_PK::MechanicsFracturedMatrix_PK(Teuchos::ParameterList& pk_tree,
                                                         const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                                         const Teuchos::RCP<State>& S,
                                                         const Teuchos::RCP<TreeVector>& soln)
  : MechanicsSmallStrain_PK(pk_tree, glist, S, soln)
{
  S_ = S;

  std::string pk_name = pk_tree.name();
  auto found = pk_name.rfind("->");
  if (found != std::string::npos) pk_name.erase(0, found + 2);

  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  ec_list_ = Teuchos::sublist(pk_list, pk_name, true);

  // We also need iscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_solver_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(ec_list_, "time integrator", true);

  // domain and primary evaluators
  domain_ = ec_list_->get<std::string>("domain name", "domain");
  displacement_key_ = Keys::getKey(domain_, "displacement");
  AddDefaultPrimaryEvaluator(S_, displacement_key_);

  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = ec_list_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("MechanicsMatrixFracture", vlist));
}


/* ******************************************************************
* Contribution of this PK for DAG.
****************************************************************** */
void
MechanicsFracturedMatrix_PK::Setup()
{
  mesh_ = S_->GetMesh("domain");
  mesh_fracture_ = S_->GetMesh("fracture");

  int d = mesh_->getSpaceDimension();

  // displacement field
  std::vector<WhetStone::SchemaItem> items;
  items.push_back(std::make_tuple(AmanziMesh::Entity_kind::NODE, WhetStone::DOF_Type::POINT, d));
  items.push_back(std::make_tuple(AmanziMesh::Entity_kind::FACE, WhetStone::DOF_Type::NORMAL_COMPONENT, 1));

  auto cvs = Operators::CreateFracturedMatrixCVS(mesh_, mesh_fracture_, items);

  if (!S_->HasRecord(displacement_key_)) {
    *S_->Require<CV_t, CVS_t>(displacement_key_, Tags::DEFAULT)
      .SetMesh(mesh_)->SetGhosted(true) = *cvs;

    eval_ = Teuchos::rcp_static_cast<EvaluatorPrimary<CV_t, CVS_t>>(
      S_->GetEvaluatorPtr(displacement_key_, Tags::DEFAULT));
  }

  MechanicsSmallStrain_PK::Setup();
}

} // namespace Mechanics
} // namespace Amanzi

