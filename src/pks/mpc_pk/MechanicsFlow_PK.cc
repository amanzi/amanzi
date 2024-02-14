/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  MPC PK

  Weak coupling of mechanics and flow PKs.
*/

#include <string>

#include "PK_Utils.hh"
#include "Transport_PK.hh"

#include "MechanicsFlow_PK.hh"

namespace Amanzi {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

/* ******************************************************************
* Constructor
****************************************************************** */
MechanicsFlow_PK::MechanicsFlow_PK(Teuchos::ParameterList& pk_tree,
                                   const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                   const Teuchos::RCP<State>& S,
                                   const Teuchos::RCP<TreeVector>& soln)
  : PK_MPCWeak(pk_tree, glist, S, soln), glist_(glist)
{
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = my_list_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("MechanicsFlow", vlist));
}


/* ******************************************************************
* Setup of PK
****************************************************************** */
void
MechanicsFlow_PK::Setup()
{
  std::string passwd("");
  hydrostatic_stress_key_ = "hydrostatic_stress";
  vol_strain_key_ = "volumetric_strain";

  S_->Require<CV_t, CVS_t>(hydrostatic_stress_key_, Tags::DEFAULT, passwd)
    .SetMesh(S_->GetMesh())
    ->SetGhosted(true)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->Require<CV_t, CVS_t>(vol_strain_key_, Tags::DEFAULT, passwd)
    .SetMesh(S_->GetMesh())
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  auto pks =
    glist_->sublist("PKs").sublist(name_).get<Teuchos::Array<std::string>>("PKs order").toVector();
  glist_->sublist("PKs")
    .sublist(pks[1])
    .sublist("physical models and assumptions")
    .set<bool>("porosity strain model", true);

  glist_->sublist("PKs")
    .sublist(pks[0])
    .sublist("physical models and assumptions")
    .set<bool>("use biot model", true);

  PK_MPCWeak::Setup();
}


/* ******************************************************************
* Extended treatment of time step in transport PK.
****************************************************************** */
bool
MechanicsFlow_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  auto pk0 = sub_pks_[0];
  bool fail = pk0->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;
  pk0->CommitStep(t_old, t_new, Tags::DEFAULT);

  auto pk1 = sub_pks_[1];
  fail = pk1->AdvanceStep(t_old, t_new, reinit);
  return fail;
}

} // namespace Amanzi
