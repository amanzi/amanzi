/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "energy_bc_factory.hh"
#include "advection_factory.hh"

#include "composite_vector_function.hh"
#include "composite_vector_function_factory.hh"

#include "bdf1_time_integrator.hh"

#include "iem.hh"
#include "eos.hh"

#include "eos_evaluator.hh"
#include "iem_evaluator.hh"
#include "two_phase_energy_evaluator.hh"
#include "enthalpy_evaluator.hh"
#include "thermal_conductivity_twophase_evaluator.hh"
#include "two_phase.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory<TwoPhase> TwoPhase::reg_("two-phase energy");


// -------------------------------------------------------------
// Setup
// -------------------------------------------------------------
void TwoPhase::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBDFBase::setup(S);
  SetupEnergy_(S);
  SetupPhysicalEvaluators_(S);
};


void TwoPhase::SetupEnergy_(const Teuchos::Ptr<State>& S) {

  // Require fields and evaluators for those fields.
  // primary variable: temperature on both cells and faces, ghosted, with 1 dof
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1); // = [1, 1]
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  S->RequireField(key_, name_)->SetMesh(S->GetMesh())
    ->SetGhosted()->SetComponents(names2, locations2, num_dofs2);

  // Get data for non-field quanitites.
  S->RequireFieldEvaluator("cell_volume");
  S->RequireScalar("density_rock");

  S->RequireField("darcy_flux")->SetMesh(S->GetMesh())->SetGhosted()
                                ->AddComponent("face", AmanziMesh::FACE, 1);

  // boundary conditions
  Teuchos::ParameterList bc_plist = plist_.sublist("boundary conditions", true);
  EnergyBCFactory bc_factory(S->GetMesh(), bc_plist);
  bc_temperature_ = bc_factory.CreateTemperature();
  bc_flux_ = bc_factory.CreateEnthalpyFlux();

  // operator for advection terms
  Operators::AdvectionFactory advection_factory;
  Teuchos::ParameterList advect_plist = plist_.sublist("Advection");
  advection_ = advection_factory.create(advect_plist, S->GetMesh());
  advection_->set_num_dofs(1);

  // operator for the diffusion terms
  bool symmetric = true;
  Teuchos::ParameterList mfd_plist = plist_.sublist("Diffusion");
  matrix_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_plist, S->GetMesh()));
  matrix_->SetSymmetryProperty(symmetric);
  matrix_->SymbolicAssembleGlobalMatrices();
  matrix_->CreateMFDmassMatrices(Teuchos::null);

  // preconditioner
  Teuchos::ParameterList mfd_pc_plist = plist_.sublist("Diffusion PC");
  preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_pc_plist, S->GetMesh()));
  preconditioner_->SetSymmetryProperty(symmetric);
  preconditioner_->SymbolicAssembleGlobalMatrices();
  preconditioner_->CreateMFDmassMatrices(Teuchos::null);
  preconditioner_->InitPreconditioner(mfd_pc_plist);
  assemble_preconditioner_ = plist_.get<bool>("assemble preconditioner", true);

  // constraint on max delta T, which kicks us out of bad iterates faster?
  dT_max_ = plist_.get<double>("maximum temperature change", 10.);
};


// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void TwoPhase::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  // Get data and evaluators needed by the PK
  // -- energy, the conserved quantity
  S->RequireField(name_)->SetMesh(S->GetMesh())->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList ee_plist = plist_.sublist("energy evaluator");
  ee_plist.set("energy key", name_);
  Teuchos::RCP<TwoPhaseEnergyEvaluator> ee =
    Teuchos::rcp(new TwoPhaseEnergyEvaluator(ee_plist));
  S->SetFieldEvaluator(name_, ee);

  // -- advection of enthalpy
  S->RequireField("enthalpy_liquid")->SetMesh(S->GetMesh())
    ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList enth_plist = plist_.sublist("enthalpy evaluator");
  enth_plist.set("enthalpy key", "enthalpy_liquid");
  Teuchos::RCP<EnthalpyEvaluator> enth =
    Teuchos::rcp(new EnthalpyEvaluator(enth_plist));
  S->SetFieldEvaluator("enthalpy_liquid", enth);

  // -- thermal conductivity
  S->RequireField("thermal_conductivity")->SetMesh(S->GetMesh())
    ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList tcm_plist =
    plist_.sublist("thermal conductivity evaluator");
  Teuchos::RCP<EnergyRelations::ThermalConductivityTwoPhaseEvaluator> tcm =
    Teuchos::rcp(new EnergyRelations::ThermalConductivityTwoPhaseEvaluator(tcm_plist));
  S->SetFieldEvaluator("thermal_conductivity", tcm);
}


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void TwoPhase::initialize(const Teuchos::Ptr<State>& S) {
  // initialize BDF stuff and physical domain stuff
  PKPhysicalBDFBase::initialize(S);

  // initialize boundary conditions
  int nfaces = S->GetMesh()->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::MFD_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  // For the boundary conditions, we currently hack in the enthalpy to
  // the boundary faces to correctly advect in a constant temperature
  // BC.  This requires density and internal energy, which in turn
  // require a model based on p,T.
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
};


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void TwoPhase::commit_state(double dt, const Teuchos::RCP<State>& S) {};


// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void TwoPhase::UpdateBoundaryConditions_() {
  for (int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::MFD_BC_NULL;
    bc_values_[n] = 0.0;
  }

  Functions::BoundaryFunction::Iterator bc;
  for (bc=bc_temperature_->begin(); bc!=bc_temperature_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
    bc_values_[f] = bc->second;
  }

  for (bc=bc_flux_->begin(); bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MFD_BC_FLUX;
    bc_values_[f] = bc->second;
  }
};


// -----------------------------------------------------------------------------
// Add a boundary marker to owned faces.
// -----------------------------------------------------------------------------
void TwoPhase::ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& temperature) {
  int nfaces = temperature->size("face");
  for (int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] == Operators::MFD_BC_DIRICHLET) {
      (*temperature)("face",f) = bc_values_[f];
    }
  }
};

bool TwoPhase::is_admissible(Teuchos::RCP<const TreeVector> up) {
  Teuchos::OSTab tab = getOSTab();
  // For some reason, wandering PKs break most frequently with an unreasonable
  // temperature.  This simply tries to catch that before it happens.
  Teuchos::RCP<const CompositeVector> temp = up->data();

  const Epetra_MultiVector& temp_v = *temp->ViewComponent("cell",false);
  double minT(0.), maxT(0.);
  int ierr = temp_v.MinValue(&minT);
  ierr |= temp_v.MaxValue(&maxT);

  if(out_.get() && includesVerbLevel(verbosity_,Teuchos::VERB_EXTREME,true)) {
    *out_ << "Admissible T? (min/max): " << minT << ",  " << maxT << std::endl;
  }

  if (ierr || minT < 200.0 || maxT > 300.0) {
    if(out_.get() && includesVerbLevel(verbosity_,Teuchos::VERB_HIGH,true)) {
      *out_ << "Energy PK is inadmissible, as it is not within bounds of constitutive models: min(T) = " << minT << ", max(T) = " << maxT << std::endl;
    }
    return false;
  }
  return true;
}

} // namespace Energy
} // namespace Amanzi
