/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "primary_variable_field_evaluator.hh"

#include "energy_bc_factory.hh"
#include "advection_factory.hh"

#include "composite_vector_function.hh"
#include "composite_vector_function_factory.hh"

#include "bdf1_time_integrator.hh"

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
// Constructor
// -------------------------------------------------------------
TwoPhase::TwoPhase(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& solution) :
    energy_plist_(plist) {
  solution_ = solution;
  SetupEnergy_(S);
  SetupPhysicalEvaluators_(S);
};


void TwoPhase::SetupEnergy_(const Teuchos::RCP<State>& S) {

  // Require fields and evaluators for those fields.
  // primary variable: temperature on both cells and faces, ghosted, with 1 dof
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1); // = [1, 1]
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  S->RequireField("temperature", "energy")->SetMesh(S->GetMesh())
    ->SetGhosted()->SetComponents(names2, locations2, num_dofs2);
  Teuchos::RCP<PrimaryVariableFieldEvaluator> temp_evaluator =
    Teuchos::rcp(new PrimaryVariableFieldEvaluator("temperature"));
  S->SetFieldEvaluator("temperature", temp_evaluator);

  // Get data for non-field quanitites.
  S->RequireFieldEvaluator("cell_volume");
  S->RequireScalar("density_rock");

  S->RequireField("darcy_flux")->SetMesh(S->GetMesh())->SetGhosted()
                                ->AddComponent("face", AmanziMesh::FACE, 1);

  // boundary conditions
  Teuchos::ParameterList bc_plist = energy_plist_.sublist("boundary conditions", true);
  EnergyBCFactory bc_factory(S->GetMesh(), bc_plist);
  bc_temperature_ = bc_factory.CreateTemperature();
  bc_flux_ = bc_factory.CreateEnthalpyFlux();

  // operator for advection terms
  Operators::AdvectionFactory advection_factory;
  Teuchos::ParameterList advect_plist = energy_plist_.sublist("Advection");
  advection_ = advection_factory.create(advect_plist, S->GetMesh());
  advection_->set_num_dofs(1);

  // operator for the diffusion terms
  bool symmetric = true;
  Teuchos::ParameterList mfd_plist = energy_plist_.sublist("Diffusion");
  matrix_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_plist, S->GetMesh()));
  matrix_->SetSymmetryProperty(symmetric);
  matrix_->SymbolicAssembleGlobalMatrices();
  matrix_->CreateMFDmassMatrices(Teuchos::null);

  // preconditioner
  Teuchos::ParameterList mfd_pc_plist = energy_plist_.sublist("Diffusion PC");
  preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_pc_plist, S->GetMesh()));
  preconditioner_->SetSymmetryProperty(symmetric);
  preconditioner_->SymbolicAssembleGlobalMatrices();
  preconditioner_->CreateMFDmassMatrices(Teuchos::null);
  preconditioner_->InitPreconditioner(mfd_pc_plist);
};


// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void TwoPhase::SetupPhysicalEvaluators_(const Teuchos::RCP<State>& S) {
  // Get data and evaluators needed by the PK
  // -- energy, the conserved quantity
  S->RequireField("energy")->SetMesh(S->GetMesh())->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList ee_plist = energy_plist_.sublist("energy evaluator");
  ee_plist.set("energy key", "energy");
  Teuchos::RCP<TwoPhaseEnergyEvaluator> ee =
    Teuchos::rcp(new TwoPhaseEnergyEvaluator(ee_plist));
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
  Teuchos::RCP<EnergyRelations::ThermalConductivityTwoPhaseEvaluator> tcm =
    Teuchos::rcp(new EnergyRelations::ThermalConductivityTwoPhaseEvaluator(tcm_plist));
  S->SetFieldEvaluator("thermal_conductivity", tcm);
}


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void TwoPhase::initialize(const Teuchos::RCP<State>& S) {
  // initial timestep size
  dt_ = energy_plist_.get<double>("initial time step", 1.);

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

  // initial conditions
  // -- Get the IC function plist.
  if (!energy_plist_.isSublist("initial condition")) {
    std::stringstream messagestream;
    messagestream << "Two Phase Energy PK has no initial condition parameter list.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  // -- Calculate the IC.
  Teuchos::ParameterList ic_plist = energy_plist_.sublist("initial condition");
  Teuchos::RCP<Field> temp_field = S->GetField("temperature", "energy");
  temp_field->Initialize(ic_plist);

  // update face temperatures as a hint?
  Teuchos::RCP<CompositeVector> temp = S->GetFieldData("temperature", "energy");
  DeriveFaceValuesFromCellValues_(S, temp);

  // initialize the timesteppper
  solution_->set_data(temp);
  atol_ = energy_plist_.get<double>("absolute error tolerance",1.0);
  rtol_ = energy_plist_.get<double>("relative error tolerance",1.0);

  if (!energy_plist_.get<bool>("strongly coupled PK", false)) {
    // -- instantiate time stepper
    Teuchos::RCP<Teuchos::ParameterList> bdf1_plist_p =
      Teuchos::rcp(new Teuchos::ParameterList(energy_plist_.sublist("time integrator")));
    time_stepper_ = Teuchos::rcp(new BDF1TimeIntegrator(this, bdf1_plist_p, solution_));
    time_step_reduction_factor_ = bdf1_plist_p->get<double>("time step reduction factor");

    // -- initialize time derivative
    Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp(new TreeVector(*solution_));
    solution_dot->PutScalar(0.0);

    // -- set initial state
    time_stepper_->set_initial_state(S->time(), solution_, solution_dot);
  }
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
// Transfer operators -- ONLY COPIES POINTERS
// -----------------------------------------------------------------------------
void TwoPhase::state_to_solution(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& solution) {
  solution->set_data(S->GetFieldData("temperature", "energy"));
};


// -----------------------------------------------------------------------------
// Transfer operators -- ONLY COPIES POINTERS
// -----------------------------------------------------------------------------
void TwoPhase::solution_to_state(const Teuchos::RCP<TreeVector>& solution,
        const Teuchos::RCP<State>& S) {
  S->SetData("temperature", "energy", solution->data());
  Teuchos::RCP<FieldEvaluator> fm = S->GetFieldEvaluator("temperature");
  Teuchos::RCP<PrimaryVariableFieldEvaluator> pri_fm =
      Teuchos::rcp_static_cast<PrimaryVariableFieldEvaluator>(fm);
  pri_fm->SetFieldAsChanged();
};


// -----------------------------------------------------------------------------
// Advance from state S to state S_next at time S.time + dt.
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// Advance from state S to state S_next at time S.time + dt.
// -----------------------------------------------------------------------------
bool TwoPhase::advance(double dt) {
  state_to_solution(S_next_, solution_);

  // take a bdf timestep
  double h = dt;
  double dt_solver;
  try {
    dt_solver = time_stepper_->time_step(h, solution_);
  } catch (Exceptions::Amanzi_exception &error) {
    if (S_next_->GetMesh()->get_comm()->MyPID() == 0) {
      std::cout << "Timestepper called error: " << error.what() << std::endl;
    }
    if (error.what() == std::string("BDF time step failed") ||
        error.what() == std::string("Cut time step")) {
      // try cutting the timestep
      dt_ = h*time_step_reduction_factor_;
      return true;
    } else {
      throw error;
    }
  }

  // commit the step as successful
  time_stepper_->commit_solution(h, solution_);
  solution_to_state(solution_, S_next_);
  commit_state(h, S_next_);

  // update the timestep size
  if (dt_solver < dt_ && dt_solver >= h) {
    // We took a smaller step than we recommended, and it worked fine (not
    // suprisingly).  Likely this was due to constraints from other PKs or
    // vis.  Do not reduce our recommendation.
  } else {
    dt_ = dt_solver;
  }

  return false;
};


// -----------------------------------------------------------------------------
// Interpolate pressure ICs on cells to ICs for lambda (faces).
// -----------------------------------------------------------------------------
void TwoPhase::DeriveFaceValuesFromCellValues_(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<CompositeVector>& temp) {
  AmanziMesh::Entity_ID_List cells;

  int f_owned = temp->size("face");
  for (int f=0; f!=f_owned; ++f) {
    cells.clear();
    S->GetMesh()->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    double face_value = 0.0;
    for (int n=0; n!=ncells; ++n) {
      face_value += (*temp)("cell",cells[n]);
    }
    (*temp)("face",f) = face_value / ncells;
  }
};


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

} // namespace Energy
} // namespace Amanzi
