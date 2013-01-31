/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */
#include "primary_variable_field_evaluator.hh"
#include "three_phase_energy_evaluator.hh"
#include "enthalpy_evaluator.hh"
#include "thermal_conductivity_threephase_evaluator.hh"
#include "three_phase.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory<ThreePhase> ThreePhase::reg_("three-phase energy");

// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void ThreePhase::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  // Get data and evaluators needed by the PK
  // -- energy, the conserved quantity
  S->RequireField(energy_key_)->SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList ee_plist = plist_.sublist("energy evaluator");
  ee_plist.set("energy key", energy_key_);
  Teuchos::RCP<ThreePhaseEnergyEvaluator> ee =
    Teuchos::rcp(new ThreePhaseEnergyEvaluator(ee_plist));
  S->SetFieldEvaluator(energy_key_, ee);

  // -- advection of enthalpy
  S->RequireField(enthalpy_key_)->SetMesh(mesh_)
    ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList enth_plist = plist_.sublist("enthalpy evaluator");
  enth_plist.set("enthalpy key", enthalpy_key_);
  Teuchos::RCP<EnthalpyEvaluator> enth =
    Teuchos::rcp(new EnthalpyEvaluator(enth_plist));
  S->SetFieldEvaluator(enthalpy_key_, enth);

  // -- thermal conductivity
  S->RequireField(conductivity_key_)->SetMesh(mesh_)
    ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList tcm_plist =
    plist_.sublist("thermal conductivity evaluator");
  Teuchos::RCP<EnergyRelations::ThermalConductivityThreePhaseEvaluator> tcm =
    Teuchos::rcp(new EnergyRelations::ThermalConductivityThreePhaseEvaluator(tcm_plist));
  S->SetFieldEvaluator(conductivity_key_, tcm);

  // require a density for rock to get total internal energy
  S->RequireScalar("density_rock");
}

} // namespace Energy
} // namespace Amanzi
