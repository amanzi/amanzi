/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code. 
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
------------------------------------------------------------------------- */


#include "bdf1_time_integrator.hh"
#include "flow_bc_factory.hh"

#include "upwind_cell_centered.hh"
#include "upwind_arithmetic_mean.hh"
#include "upwind_total_flux.hh"
#include "upwind_gravity_flux.hh"

#include "composite_vector_function.hh"
#include "composite_vector_function_factory.hh"

#include "primary_variable_field_evaluator.hh"
#include "wrm_permafrost_evaluator.hh"
#include "rel_perm_evaluator.hh"
#include "permafrost_water_content.hh"

#include "permafrost.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory<Permafrost> Permafrost::reg_("permafrost flow");


// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void Permafrost::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  // -- Absolute permeability.
  //       For now, we assume scalar permeability.  This will change.
  S->RequireField("permeability")->SetMesh(S->GetMesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("permeability");

  // -- water content, and evaluator
  S->RequireField("water_content")->SetMesh(S->GetMesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList wc_plist = plist_.sublist("water content evaluator");
  Teuchos::RCP<PermafrostWaterContent> wc =
      Teuchos::rcp(new PermafrostWaterContent(wc_plist));
  S->SetFieldEvaluator("water_content", wc);

  // -- Water retention evaluators, for saturation and rel perm.
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::BOUNDARY_FACE;
  names2[0] = "cell";
  names2[1] = "boundary_face";

  // -- rel perm on cells + boundary faces
  S->RequireField("relative_permeability")->SetMesh(S->GetMesh())->SetGhosted()
      ->AddComponents(names2,locations2,num_dofs2);

  // -- This setup is a little funky -- we use four evaluators to capture the physics.
  Teuchos::ParameterList wrm_plist = plist_.sublist("water retention evaluator");
  Teuchos::RCP<FlowRelations::WRMPermafrostEvaluator> wrm =
      Teuchos::rcp(new FlowRelations::WRMPermafrostEvaluator(wrm_plist));
  S->SetFieldEvaluator("saturation_liquid", wrm);
  S->SetFieldEvaluator("saturation_gas", wrm);
  S->SetFieldEvaluator("saturation_ice", wrm);

  // -- the rel perm evaluator, also with the same underlying WRM.
  Teuchos::RCP<FlowRelations::RelPermEvaluator> rel_perm_evaluator =
      Teuchos::rcp(new FlowRelations::RelPermEvaluator(wrm_plist, wrm->get_WRMs()));
  S->SetFieldEvaluator("relative_permeability", rel_perm_evaluator);

  // -- Liquid density and viscosity for the transmissivity.
  S->RequireField("molar_density_liquid")->SetMesh(S->GetMesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("molar_density_liquid");

  S->RequireField("viscosity_liquid")->SetMesh(S->GetMesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("viscosity_liquid");

  // -- liquid mass density for the gravity fluxes
  S->RequireField("mass_density_liquid")->SetMesh(S->GetMesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("mass_density_liquid"); // simply picks up the molar density one.
}

} // namespace
} // namespace
