/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Process kernel for energy equation for Richard's flow.
------------------------------------------------------------------------- */

#include "eos_evaluator.hh"
#include "iem_evaluator.hh"
#include "thermal_conductivity_twophase_evaluator.hh"
#include "two_phase_energy_evaluator.hh"
#include "enthalpy_evaluator.hh"
#include "energy_bc_factory.hh"

#include "two_phase.hh"

namespace Amanzi {
namespace Energy {

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
TwoPhase::TwoPhase(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<TreeVector>& solution) :
    PKDefaultBase(plist, solution),
    EnergyBase(plist, solution) {
  if (!plist_.isParameter("flux key")) plist_.set("flux key", "darcy_flux");
}

// -------------------------------------------------------------
// Create the physical evaluators for energy, enthalpy, thermal
// conductivity, and any sources.
// -------------------------------------------------------------
void TwoPhase::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  // Get data and evaluators needed by the PK
  // -- energy, the conserved quantity
  S->RequireField(energy_key_)->SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList ee_plist = plist_.sublist("energy evaluator");
  ee_plist.set("energy key", energy_key_);
  Teuchos::RCP<TwoPhaseEnergyEvaluator> ee =
    Teuchos::rcp(new TwoPhaseEnergyEvaluator(ee_plist));
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
  Teuchos::RCP<EnergyRelations::ThermalConductivityTwoPhaseEvaluator> tcm =
    Teuchos::rcp(new EnergyRelations::ThermalConductivityTwoPhaseEvaluator(tcm_plist));
  S->SetFieldEvaluator(conductivity_key_, tcm);

  // require a density for rock to get total internal energy
  S->RequireScalar("density_rock");
}


// -------------------------------------------------------------
// Initialize the needed models to plug in enthalpy.
// -------------------------------------------------------------
void TwoPhase::initialize(const Teuchos::Ptr<State>& S) {
  // Call the base class's initialize.
  EnergyBase::initialize(S);

  // For the boundary conditions, we currently hack in the enthalpy to
  // the boundary faces to correctly advect in a Dirichlet temperature
  // BC.  This requires density and internal energy, which in turn
  // require a model based on p,T.
  // This will be removed once boundary faces are implemented.
  Teuchos::RCP<FieldEvaluator> eos_fe = S->GetFieldEvaluator("molar_density_liquid");
  Teuchos::RCP<Relations::EOSEvaluator> eos_eval =
    Teuchos::rcp_dynamic_cast<Relations::EOSEvaluator>(eos_fe);
  ASSERT(eos_eval != Teuchos::null);
  eos_liquid_ = eos_eval->get_EOS();

  Teuchos::RCP<FieldEvaluator> iem_fe = S->GetFieldEvaluator("internal_energy_liquid");
  Teuchos::RCP<EnergyRelations::IEMEvaluator> iem_eval =
    Teuchos::rcp_dynamic_cast<EnergyRelations::IEMEvaluator>(iem_fe);
  ASSERT(iem_eval != Teuchos::null);
  iem_liquid_ = iem_eval->get_IEM();

}


// -------------------------------------------------------------
// Plug enthalpy into the boundary faces manually.
// This will be removed once boundary faces exist.
// -------------------------------------------------------------
void TwoPhase::ApplyDirichletBCsToEnthalpy_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& enth) {

  // put the boundary fluxes in faces for Dirichlet BCs.
  // NOTE this boundary flux is in enthalpy, and
  // h = n(T,p) * u_l(T) + p_l
  const Epetra_MultiVector& pres = *S->GetFieldData("pressure")
      ->ViewComponent("face",false);
  const Epetra_MultiVector& flux = *S->GetFieldData(flux_key_)
      ->ViewComponent("face",false);

  Epetra_MultiVector& enth_f = *enth->ViewComponent("face",false);

  for (Functions::BoundaryFunction::Iterator bc = bc_temperature_->begin();
       bc!=bc_temperature_->end(); ++bc) {
    int f = bc->first;
    double p = pres[0][f];
    double T = bc->second;
    double dens = eos_liquid_->MolarDensity(T,p);
    double int_energy = iem_liquid_->InternalEnergy(T);
    double enthalpy = int_energy + p/dens;

    enth_f[0][f] = enthalpy * fabs(flux[0][f]);
  }

}



} // namespace
} // namespace
