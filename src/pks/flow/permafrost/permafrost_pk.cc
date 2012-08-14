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

#include "primary_variable_field_model.hh"
#include "wrm_permafrost_model.hh"
#include "wrm_standard_model.hh"
#include "wrm_ice_water_model.hh"
#include "rel_perm_model.hh"
#include "permafrost_water_content.hh"

#include "permafrost.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory<Permafrost> Permafrost::reg_("permafrost flow");


// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
Permafrost::Permafrost(Teuchos::ParameterList& flow_plist,
                       const Teuchos::RCP<State>& S,
                       const Teuchos::RCP<TreeVector>& solution) {
  flow_plist_ = flow_plist;
  solution_ = solution;

  // Creation is done in two parts.
  // -- Pieces common to all flow models.
  SetupRichardsFlow_(S);

  // -- Pieces specific to the permafrost model.
  SetupPhysicalModels_(S);
};

// -------------------------------------------------------------
// Create the physical models for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void Permafrost::SetupPhysicalModels_(const Teuchos::RCP<State>& S) {
  // -- Absolute permeability.
  //       For now, we assume scalar permeability.  This will change.
  S->RequireField("permeability")->SetMesh(S->Mesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldModel("permeability");

  // -- water content, and model
  S->RequireField("water_content")->SetMesh(S->Mesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList wc_plist = flow_plist_.sublist("water content model");
  Teuchos::RCP<PermafrostWaterContent> wc =
      Teuchos::rcp(new PermafrostWaterContent(wc_plist));
  S->SetFieldModel("water_content", wc);

  // -- Water retention models, for saturation and rel perm.
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  S->RequireField("relative_permeability")->SetMesh(S->Mesh())->SetGhosted()
                    ->AddComponents(names2, locations2, num_dofs2);

  // -- This setup is a little funky -- we use four models to capture the physics.

  //    This follows the permafrost notes... the first model,
  //    WRMStandardModel, is the usual S(p_atm - p_liquid), which in the
  //    permafrost model is given by 1/B.  The second model, WRMIceWaterModel,
  //    is given by S( gamma * T' ), or 1/A.  The third model depends upon
  //    these two models, and evaluates s_g, s_l, and s_i from 1/A and 1/B.
  //
  //    The fourth model provides Krel(s_liquid), which is NOT the same as
  //    Krel(p_atm - p_liquid).

  // Model 3.
  Teuchos::ParameterList wrm_plist = flow_plist_.sublist("water retention model");
  Teuchos::RCP<FlowRelations::WRMPermafrostModel> wrm =
      Teuchos::rcp(new FlowRelations::WRMPermafrostModel(wrm_plist));
  S->SetFieldModel("saturation_liquid", wrm);
  S->SetFieldModel("saturation_gas", wrm);
  S->SetFieldModel("saturation_ice", wrm);

  // Model 1.
  Teuchos::ParameterList Aplist;
  std::string Akey = wrm_plist.get<string>("1/A key", "wrm_permafrost_one_on_A");
  Aplist.set("saturation key", Akey);
  Aplist.set("calculate minor saturation", false);
  ASSERT(wrm_plist.isSublist("WRM parameters"));
  Aplist.set("WRM parameters", wrm_plist.sublist("WRM parameters"));
  Teuchos::RCP<FlowRelations::WRMStandardModel> wrm_A =
      Teuchos::rcp(new FlowRelations::WRMStandardModel(Aplist));
  S->SetFieldModel(Akey, wrm_A);

  // Model 2.  Constructed using the same underlying vanGenuchten model as model 1.
  Teuchos::ParameterList Bplist;
  std::string Bkey = wrm_plist.get<string>("1/B key", "wrm_permafrost_one_on_B");
  Bplist.set("saturation key", Bkey);
  Bplist.set("calculate minor saturation", false);
  Teuchos::RCP<FlowRelations::WRMIceLiquidModel> wrm_B =
      Teuchos::rcp(new FlowRelations::WRMIceLiquidModel(Bplist, wrm_A->get_WRM()));
  S->SetFieldModel(Bkey, wrm_B);

  // Model 4, the rel perm model, also with the same underlying model.
  Teuchos::RCP<FlowRelations::RelPermModel> rel_perm_model =
      Teuchos::rcp(new FlowRelations::RelPermModel(wrm_plist, wrm_A->get_WRM()));
  S->SetFieldModel("relative_permeability", rel_perm_model);


  // -- Liquid density and viscosity for the transmissivity.
  S->RequireField("molar_density_liquid")->SetMesh(S->Mesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldModel("molar_density_liquid");

  S->RequireField("viscosity_liquid")->SetMesh(S->Mesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldModel("viscosity_liquid");
}

} // namespace
} // namespace
