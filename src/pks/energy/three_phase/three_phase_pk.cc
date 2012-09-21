/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "three_phase_energy_evaluator.hh"
#include "enthalpy_evaluator.hh"
#include "thermal_conductivity_threephase_evaluator.hh"
#include "three_phase.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory<ThreePhase> ThreePhase::reg_("three-phase energy");


// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
ThreePhase::ThreePhase(Teuchos::ParameterList& plist,
                       const Teuchos::RCP<State>& S,
                       const Teuchos::RCP<TreeVector>& solution) {
  energy_plist_ = plist;
  solution_ = solution;
  SetupEnergy_(S);
  SetupPhysicalEvaluators_(S);
};


// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void ThreePhase::SetupPhysicalEvaluators_(const Teuchos::RCP<State>& S) {
  // Get data and evaluators needed by the PK
  // -- energy, the conserved quantity
  S->RequireField("energy")->SetMesh(S->GetMesh())->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList ee_plist = energy_plist_.sublist("energy evaluator");
  ee_plist.set("energy key", "energy");
  Teuchos::RCP<ThreePhaseEnergyEvaluator> ee =
    Teuchos::rcp(new ThreePhaseEnergyEvaluator(ee_plist));
  S->SetFieldEvaluator("energy", ee);

  // -- advection of enthalpy
  S->RequireField("enthalpy_liquid")->SetMesh(S->GetMesh())
    ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList enth_plist = energy_plist_.sublist("enthalpy evaluator");
  enth_plist.set("enthalpy key", "enthalpy_liquid");
  Teuchos::RCP<EnthalpyEvaluator> enth =
    Teuchos::rcp(new EnthalpyEvaluator(enth_plist));
  S->SetFieldEvaluator("enthalpy_liquid", enth);

  // -- thermal conductivity
  S->RequireField("thermal_conductivity")->SetMesh(S->GetMesh())
    ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList tcm_plist =
    energy_plist_.sublist("thermal conductivity evaluator");
  Teuchos::RCP<EnergyRelations::ThermalConductivityThreePhaseEvaluator> tcm =
    Teuchos::rcp(new EnergyRelations::ThermalConductivityThreePhaseEvaluator(tcm_plist));
  S->SetFieldEvaluator("thermal_conductivity", tcm);
}

} // namespace Energy
} // namespace Amanzi
