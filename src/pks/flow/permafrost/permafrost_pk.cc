/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code. 
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
------------------------------------------------------------------------- */


#include "flow_bc_factory.hh"

#include "upwind_cell_centered.hh"
#include "upwind_arithmetic_mean.hh"
#include "upwind_total_flux.hh"
#include "upwind_gravity_flux.hh"

#include "CompositeVectorFunction.hh"
#include "CompositeVectorFunctionFactory.hh"

#include "primary_variable_field_evaluator.hh"
#include "wrm_permafrost_evaluator.hh"
#include "rel_perm_evaluator.hh"

#include "permafrost.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void Permafrost::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  // -- Absolute permeability.
  //       For now, we assume scalar permeability.  This will change.
  S->RequireField(getKey(domain_,"permeability"))->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(getKey(domain_,"permeability"));
  
  S->RequireField(conserved_key_)->SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(conserved_key_);

  // -- Water retention evaluators, for saturation and rel perm.
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::BOUNDARY_FACE;
  names2[0] = "cell";
  names2[1] = "boundary_face";

  // -- rel perm on cells + boundary faces

  S->RequireField(coef_key_)->SetMesh(mesh_)->SetGhosted()
      ->AddComponents(names2,locations2,num_dofs2);
 
  // -- This setup is a little funky -- we use four evaluators to capture the physics.
  Teuchos::ParameterList wrm_plist = plist_->sublist("water retention evaluator");
  Teuchos::RCP<FlowRelations::WRMPermafrostEvaluator> wrm =
      Teuchos::rcp(new FlowRelations::WRMPermafrostEvaluator(wrm_plist));

  
  S->SetFieldEvaluator(getKey(domain_,"saturation_liquid"), wrm);

  S->SetFieldEvaluator(getKey(domain_,"saturation_gas"), wrm);

  S->SetFieldEvaluator(getKey(domain_,"saturation_ice"), wrm);

  // -- the rel perm evaluator, also with the same underlying WRM.
  wrm_plist.set<double>("permeability rescaling", perm_scale_);
  Teuchos::RCP<FlowRelations::RelPermEvaluator> rel_perm_evaluator =
      Teuchos::rcp(new FlowRelations::RelPermEvaluator(wrm_plist, wrm->get_WRMs()));
  wrms_ = wrm->get_WRMs();
  

//I-COMMENTED  S->SetFieldEvaluator("relative_permeability", rel_perm_evaluator);
 S->SetFieldEvaluator(coef_key_, rel_perm_evaluator);

  // -- Liquid density and viscosity for the transmissivity.
  S->RequireField(getKey(domain_,"molar_density_liquid"))->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(getKey(domain_,"molar_density_liquid"));

  /* S->RequireField("viscosity_liquid")->SetMesh(S->GetMesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("viscosity_liquid");
  */
  // -- liquid mass density for the gravity fluxes
  S->RequireField(getKey(domain_,"mass_density_liquid"))->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(getKey(domain_,"mass_density_liquid")); // simply picks up the molar density one.

}

} // namespace
} // namespace
