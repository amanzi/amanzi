/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include "bdf2_time_integrator.hh"
#include "flow_bc_factory.hh"
#include "eos_factory.hh"

#include "richards.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory<Richards> Richards::reg_("richards flow");

// constructor
Richards::Richards(Teuchos::ParameterList& flow_plist, const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& solution) :
    flow_plist_(flow_plist) {


  // just the extras...
  // data layouts for fields
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  // -- primary variable: pressure on both cells and faces, ghosted, with 1 dof
  std::vector< std::vector<std::string> > subfield_names(2);
  subfield_names[0].resize(1);
  subfield_names[0][0] = "pressure";
  subfield_names[1].resize(1);
  subfield_names[1][0] = "pressure_lambda";
  S->RequireField("pressure", "flow", names2, locations2, 1, true);
  Teuchos::RCP<CompositeVector> pressure = S->GetFieldData("pressure", "flow");
  pressure->set_subfield_names(subfield_names);

  // update the tree-vector solution with flow's primary variable
  solution->set_data(pressure);
  solution_ = solution;

  // -- secondary variables
  S->RequireField("darcy_flux", "flow", AmanziMesh::FACE, 1, true);
  S->RequireField("darcy_velocity", "flow", AmanziMesh::CELL, 3, false);

  S->RequireField("saturation_liquid", "flow", AmanziMesh::CELL, 1, true);
  S->RequireField("density_liquid", "flow", AmanziMesh::CELL, 1, true);
  S->RequireField("molar_density_liquid", "flow", AmanziMesh::CELL, 1, true);
  S->RequireField("viscosity_liquid", "flow", AmanziMesh::CELL, 1, true);

  S->RequireField("saturation_gas", "flow", AmanziMesh::FACE, 1, true);
  S->RequireField("density_gas", "flow", AmanziMesh::CELL, 1, true);
  S->RequireField("mol_frac_gas", "flow", AmanziMesh::CELL, 1, true);
  S->RequireField("molar_density_gas", "flow", AmanziMesh::CELL, 1, true);

  S->RequireField("permeability", "flow", AmanziMesh::CELL, 1, true);
  S->RequireField("relative_permeability", "flow", AmanziMesh::CELL, 1, true);
  S->RequireScalar("atmospheric_pressure", "flow");

  // -- independent variables not owned by this PK
  S->RequireField("porosity", AmanziMesh::CELL, 1, true);
  S->RequireField("cell_volume", AmanziMesh::CELL, 1, true);
  S->RequireField("temperature", AmanziMesh::CELL, 1, true);
  S->RequireConstantVector("gravity", S->mesh()->space_dimension());

  // -- work vectors
  S->RequireField("rel_perm_faces", "flow", AmanziMesh::FACE, 1, true);
  S->GetRecord("rel_perm_faces","flow")->set_io_vis(false);

  // abs perm tensor
  variable_abs_perm_ = false; // currently not implemented, but may eventually want a model
  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  K_.resize(c_owned);
  for (int c=0; c!=c_owned; ++c) K_[c].init(S->mesh()->space_dimension(),1);

  // constitutive relations
  Teuchos::ParameterList water_eos_plist = flow_plist_.sublist("Water EOS");
  FlowRelations::EOSFactory eos_factory;
  eos_liquid_ = eos_factory.createEOS(water_eos_plist);

  Teuchos::ParameterList eos_gas_plist = flow_plist_.sublist("Vapor and Gas EOS");
  eos_gas_ = Teuchos::rcp(new FlowRelations::EOSVaporInGas(eos_gas_plist));

  Teuchos::ParameterList wrm_plist = flow_plist_.sublist("Water Retention Model");
  wrm_ = Teuchos::rcp(new FlowRelations::WRMVanGenuchten(wrm_plist));

  // boundary conditions
  Teuchos::ParameterList bc_plist = flow_plist_.sublist("boundary conditions", true);
  FlowBCFactory bc_factory(S->mesh(), bc_plist);
  bc_pressure_ = bc_factory.CreatePressure();
  bc_flux_ = bc_factory.CreateMassFlux();
  // head boundary conditions not yet implemented, as they depend on single density
  // bc_head_ = bc_factory.CreateStaticHead(0.0, 0.0, g);

  // operator for the diffusion terms
  Teuchos::ParameterList mfd_plist = flow_plist_.sublist("Diffusion");
  matrix_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_plist, S->mesh()));
  matrix_->SetSymmetryProperty(true);
  matrix_->SymbolicAssembleGlobalMatrices();

  // Relative permeability method
  string method_name = mfd_plist.get<string>("Relative permeability method", "Upwind with gravity");
  if (method_name == "Upwind with gravity") {
    Krel_method_ = FLOW_RELATIVE_PERM_UPWIND_GRAVITY;
  } else if (method_name == "Cell centered") {
    Krel_method_ = FLOW_RELATIVE_PERM_CENTERED;
  } else if (method_name == "Upwind with Darcy flux") {
    Krel_method_ = FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX;
  } else if (method_name == "Arithmetic mean") {
    Krel_method_ = FLOW_RELATIVE_PERM_ARITHMETIC_MEAN;
  }

  // preconditioner for the NKA system
  Teuchos::ParameterList mfd_pc_plist = flow_plist_.sublist("Diffusion PC");
  preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_pc_plist, S->mesh()));
  preconditioner_->SetSymmetryProperty(true);
  preconditioner_->SymbolicAssembleGlobalMatrices();
  Teuchos::ParameterList mfd_pc_ml_plist = mfd_pc_plist.sublist("ML Parameters");
  preconditioner_->InitMLPreconditioner(mfd_pc_ml_plist);
};

// -- Initialize owned (dependent) variables.
void Richards::initialize(const Teuchos::RCP<State>& S) {

  // initial timestep size
  dt_ = flow_plist_.get<double>("Initial time step", 1.);

  // initialize pressure?
  if (flow_plist_.isParameter("Constant pressure")) {
    double p = flow_plist_.get<double>("Constant pressure");
    S->GetFieldData("pressure", "flow")->PutScalar(p);
    S->GetRecord("pressure", "flow")->set_initialized();
  }

  // update secondary variables for IC
  UpdateSecondaryVariables_(S);

  // declare secondary variables initialized
  S->GetRecord("saturation_liquid","flow")->set_initialized();
  S->GetRecord("saturation_gas","flow")->set_initialized();
  S->GetRecord("density_liquid","flow")->set_initialized();
  S->GetRecord("viscosity_liquid","flow")->set_initialized();
  S->GetRecord("density_gas","flow")->set_initialized();
  S->GetRecord("molar_density_liquid","flow")->set_initialized();
  S->GetRecord("molar_density_gas","flow")->set_initialized();
  S->GetRecord("mol_frac_gas","flow")->set_initialized();
  S->GetRecord("relative_permeability","flow")->set_initialized();
  S->GetRecord("darcy_flux", "flow")->set_initialized();
  S->GetRecord("darcy_velocity", "flow")->set_initialized();

  S->GetFieldData("rel_perm_faces","flow")->PutScalar(1.0);
  S->GetRecord("rel_perm_faces","flow")->set_initialized();

  // initialize boundary conditions
  int nfaces = S->mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::MFD_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  double time = S->time();
  bc_pressure_->Compute(time);
  //bc_head_->Compute(time);
  bc_flux_->Compute(time);
  UpdateBoundaryConditions_();

  // initialize the timesteppper
  state_to_solution(S, solution_);
  if (!flow_plist_.get<bool>("Strongly Coupled PK", false)) {
    // model evaluator params
    // -- tolerances
    atol_ = flow_plist_.get<double>("Absolute error tolerance",1e-5);
    rtol_ = flow_plist_.get<double>("Relative error tolerance",1e-5);

    // -- instantiate time stepper
    Teuchos::RCP<Teuchos::ParameterList> bdf2_plist_p =
      Teuchos::rcp(new Teuchos::ParameterList(flow_plist_.sublist("Time integrator")));
    time_stepper_ = Teuchos::rcp(new BDF2TimeIntegrator(this, bdf2_plist_p, solution_));

    // -- initialize time derivative
    Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp( new TreeVector(*solution_));
    solution_dot->PutScalar(0.0);

    // -- set initial state
    time_stepper_->set_initial_state(S->time(), solution_, solution_dot);
  }
};

// Pointer copy of state to solution
void Richards::state_to_solution(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& solution) {
  solution->set_data(S->GetFieldData("pressure", "flow"));
};

// Pointer copy concentration fields from solution vector into state.  Used within
// compute_f() of strong couplers to set the current iterate in the state for
// use with other PKs.
void Richards::solution_to_state(const Teuchos::RCP<TreeVector>& solution,
        const Teuchos::RCP<State>& S) {
  S->SetData("pressure", "flow", solution->data());
};

// -- Advance from state S0 to state S1 at time S0.time + dt.
bool Richards::advance(double dt) {
  state_to_solution(S_next_, solution_);

  // take the bdf timestep
  double h = dt;
  time_stepper_->time_step(h, solution_);
  time_stepper_->commit_solution(h, solution_);
  return false;
};

// -- commit the state... not sure how/when this gets used
void Richards::commit_state(double dt, const Teuchos::RCP<State>& S) {
  // update the flux
  Teuchos::RCP<CompositeVector> flux = S->GetFieldData("darcy_flux", "flow");
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
  DeriveDarcyFlux_(*pres, flux);
};


// -- update diagnostics -- used prior to vis
void Richards::calculate_diagnostics(const Teuchos::RCP<State>& S) {
  // update the cell velocities
  Teuchos::RCP<CompositeVector> velocity = S->GetFieldData("darcy_velocity", "flow");
  Teuchos::RCP<const CompositeVector> flux = S->GetFieldData("darcy_flux");
  DeriveDarcyVelocity_(*flux, velocity);
};


// relative permeability methods
void Richards::UpdatePermeabilityData_(const Teuchos::RCP<State>& S) {
  SetAbsolutePermeabilityTensor_(S);

  Teuchos::RCP<const CompositeVector> rel_perm = S->GetFieldData("relative_permeability");
  rel_perm->ScatterMasterToGhosted();

  Teuchos::RCP<const CompositeVector> dens = S->GetFieldData("density_liquid");
  Teuchos::RCP<const CompositeVector> visc = S->GetFieldData("viscosity_liquid");
  Teuchos::RCP<CompositeVector> rel_perm_faces = S->GetFieldData("rel_perm_faces", "flow");

  if (Krel_method_ == FLOW_RELATIVE_PERM_CENTERED) {
    // symmetric method, no faces needed
    for (int c=0; c!=K_.size(); ++c) {
      K_[c] *= (*rel_perm)(c) * (*dens)(c) / (*visc)(c);
    }
  } else {
    // faces needed
    if (Krel_method_ == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) {
      CalculateRelativePermeabilityUpwindGravity_(S, *rel_perm, rel_perm_faces);
    } else if (Krel_method_ == FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX) {
      Teuchos::RCP<const CompositeVector> flux = S_->GetFieldData("darcy_flux");
      CalculateRelativePermeabilityUpwindFlux_(S, *flux, *rel_perm,
              rel_perm_faces);
    } else if (Krel_method_ == FLOW_RELATIVE_PERM_ARITHMETIC_MEAN) {
      CalculateRelativePermeabilityArithmeticMean_(S, *rel_perm, rel_perm_faces);
    }
    // update K with just rho/mu
    for (int c=0; c!=K_.size(); ++c) {
      K_[c] *= (*dens)(c) / (*visc)(c);
    }
  }
};

/* ******************************************************************
* Defines upwinded relative permeabilities for faces using gravity.
****************************************************************** */
void Richards::CalculateRelativePermeabilityUpwindGravity_(const Teuchos::RCP<State>& S,
        const CompositeVector& rel_perm_cells,
        const Teuchos::RCP<CompositeVector>& rel_perm_faces) {
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Teuchos::RCP<const Epetra_Vector> g_vec = S->GetConstantVectorData("gravity");
  AmanziGeometry::Point gravity(g_vec->MyLength());
  for (int i=0; i!=g_vec->MyLength(); ++i) gravity[i] = (*g_vec)[i];

  rel_perm_faces->PutScalar(0.0);
  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    S->mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);

    AmanziGeometry::Point Kgravity = K_[c] * gravity;

    for (int n=0; n!=faces.size(); ++n) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = S->mesh()->face_normal(f);
      if ((normal * Kgravity) * dirs[n] >= 0.0) {
        (*rel_perm_faces)(f) = rel_perm_cells(c);
      } else if (bc_markers_[f] != Operators::MFD_BC_NULL) {
        (*rel_perm_faces)(f) = rel_perm_cells(c);
      }
    }
  }
}


/* ******************************************************************
* Defines upwinded relative permeabilities for faces using a given flux.
****************************************************************** */
void Richards::CalculateRelativePermeabilityUpwindFlux_(const Teuchos::RCP<State>& S,
        const CompositeVector& flux, const CompositeVector& rel_perm_cells,
        const Teuchos::RCP<CompositeVector>& rel_perm_faces) {
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  rel_perm_faces->PutScalar(0.0);
  int c_owned = S->mesh()->count_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c=0; c!=c_owned; ++c) {
    S->mesh()->cell_get_faces_and_dirs(c, &faces, &dirs);

    for (int n=0; n!=faces.size(); ++n) {
      int f = faces[n];
      if (flux(n) * dirs[n] >= 0.0) {
        (*rel_perm_faces)(f) = rel_perm_cells(c);
      } else if (bc_markers_[f] != Operators::MFD_BC_NULL) {
        (*rel_perm_faces)(f) = rel_perm_cells(c);
      }
    }
  }
}


/* ******************************************************************
* Defines relative permeabilities for faces via arithmetic averaging.
****************************************************************** */
void Richards::CalculateRelativePermeabilityArithmeticMean_(const Teuchos::RCP<State>& S,
        const CompositeVector& rel_perm_cells,
        const Teuchos::RCP<CompositeVector>& rel_perm_faces) {
  AmanziMesh::Entity_ID_List cells;

  rel_perm_faces->PutScalar(0.0);
  int f_owned = S->mesh()->count_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f=0; f!=f_owned; ++f) {
    S->mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
    for (int n=0; n!=cells.size(); ++n) (*rel_perm_faces)(f) += rel_perm_cells(cells[n]);
    (*rel_perm_faces)(f) /= cells.size();
  }
}


/* ******************************************************************
 * Add a boundary marker to used faces.
 ****************************************************************** */
void Richards::UpdateBoundaryConditions_() {
  for (int n=0; n<bc_markers_.size(); n++) {
    bc_markers_[n] = Operators::MFD_BC_NULL;
    bc_values_[n] = 0.0;
  }

  BoundaryFunction::Iterator bc;
  for (bc=bc_pressure_->begin(); bc!=bc_pressure_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
    bc_values_[f] = bc->second;
  }

  // for (bc=bc_head_->begin(); bc!=bc_head_->end(); ++bc) {
  //   int f = bc->first;
  //   bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
  //   bc_values_[f] = bc->second;
  // }

  for (bc=bc_flux_->begin(); bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MFD_BC_FLUX;
    bc_values_[f] = bc->second;
  }
};


/* ******************************************************************
 * Add a boundary marker to owned faces.
 ****************************************************************** */
void Richards::ApplyBoundaryConditions_(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<CompositeVector>& pres) {
  int nfaces = S->mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] == Operators::MFD_BC_DIRICHLET) {
      (*pres)("face",0,f) = bc_values_[f];
    }
  }
};


} // namespace
} // namespace
