/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code.
License: BSD
Authors: Neil Carlson (version 1)
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Ethan Coon (ATS version) (ecoon@lanl.gov)
------------------------------------------------------------------------- */
#include "boost/math/special_functions/fpclassify.hpp"

#include "bdf1_time_integrator.hh"
#include "flow_bc_factory.hh"

#include "matrix_mfd.cc"
#include "upwind_cell_centered.hh"
#include "upwind_arithmetic_mean.hh"
#include "upwind_total_flux.hh"
#include "upwind_gravity_flux.hh"

#include "composite_vector_function.hh"
#include "composite_vector_function_factory.hh"

#include "wrm_richards_evaluator.hh"
#include "rel_perm_evaluator.hh"
#include "richards_water_content.hh"

#include "richards.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory<Richards> Richards::reg_("richards flow");


// -------------------------------------------------------------
// Setup data
// -------------------------------------------------------------
void Richards::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBDFBase::setup(S);
  SetupRichardsFlow_(S);
  SetupPhysicalEvaluators_(S);
};


// -------------------------------------------------------------
// Pieces of the construction process that are common to all
// Richards-like PKs.
// -------------------------------------------------------------
void Richards::SetupRichardsFlow_(const Teuchos::Ptr<State>& S) {

  // Require fields and evaluators for those fields.
  // -- primary variable: pressure on both cells and faces, ghosted, with 1 dof
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  S->RequireField(key_, name_)->SetMesh(S->GetMesh())->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);

  // -- secondary variables, no evaluator used
  S->RequireField("darcy_flux", name_)->SetMesh(S->GetMesh())->SetGhosted()
                                ->SetComponent("face", AmanziMesh::FACE, 1);
  S->RequireField("darcy_velocity", name_)->SetMesh(S->GetMesh())->SetGhosted()
                                ->SetComponent("cell", AmanziMesh::CELL, 3);

  // Get data for non-field quanitites.
  S->RequireFieldEvaluator("cell_volume");
  S->RequireGravity();
  S->RequireScalar("atmospheric_pressure");

  // Create the absolute permeability tensor.
  int c_owned = S->GetMesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  K_ = Teuchos::rcp(new std::vector<WhetStone::Tensor>(c_owned));
  for (int c=0; c!=c_owned; ++c) {
    (*K_)[c].init(S->GetMesh()->space_dimension(),1);
  }

  // Create the boundary condition data structures.
  Teuchos::ParameterList bc_plist = plist_.sublist("boundary conditions", true);
  FlowBCFactory bc_factory(S->GetMesh(), bc_plist);
  bc_pressure_ = bc_factory.CreatePressure();
  bc_flux_ = bc_factory.CreateMassFlux();

  // Create the upwinding method
  S->RequireField("numerical_rel_perm", name_)->SetMesh(S->GetMesh())->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);
  S->GetField("numerical_rel_perm",name_)->set_io_vis(false);
  string method_name = plist_.get<string>("relative permeability method", "upwind with gravity");
  bool symmetric = false;
  if (method_name == "upwind with gravity") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindGravityFlux(name_,
            "relative_permeability", "numerical_rel_perm", K_));
    Krel_method_ = FLOW_RELATIVE_PERM_UPWIND_GRAVITY;
  } else if (method_name == "cell centered") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindCellCentered(name_,
            "relative_permeability", "numerical_rel_perm"));
    symmetric = true;
    Krel_method_ = FLOW_RELATIVE_PERM_CENTERED;
  } else if (method_name == "upwind with Darcy flux") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
            "relative_permeability", "numerical_rel_perm", "darcy_flux"));
    Krel_method_ = FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX;
  } else if (method_name == "arithmetic mean") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindArithmeticMean(name_,
            "relative_permeability", "numerical_rel_perm"));
    Krel_method_ = FLOW_RELATIVE_PERM_ARITHMETIC_MEAN;
  } else {
    std::stringstream messagestream;
    messagestream << "Richards FLow PK has no upwinding method named: " << method_name;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  // operator for the diffusion terms
  Teuchos::ParameterList mfd_plist = plist_.sublist("Diffusion");
  matrix_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_plist, S->GetMesh()));
  matrix_->SetSymmetryProperty(symmetric);
  matrix_->SymbolicAssembleGlobalMatrices();

  // preconditioner for the NKA system
  Teuchos::ParameterList mfd_pc_plist = plist_.sublist("Diffusion PC");
  preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_pc_plist, S->GetMesh()));
  preconditioner_->SetSymmetryProperty(symmetric);
  preconditioner_->SymbolicAssembleGlobalMatrices();
  preconditioner_->InitPreconditioner(mfd_pc_plist);
}

// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void Richards::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  // -- Absolute permeability.
  //       For now, we assume scalar permeability.  This will change.
  S->RequireField("permeability")->SetMesh(S->GetMesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("permeability");

  // -- water content, and evaluator
  S->RequireField("water_content")->SetMesh(S->GetMesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList wc_plist = plist_.sublist("water content evaluator");
  Teuchos::RCP<RichardsWaterContent> wc = Teuchos::rcp(new RichardsWaterContent(wc_plist));
  S->SetFieldEvaluator("water_content", wc);

  // -- Water retention evaluators, for saturation and rel perm.
  S->RequireField("relative_permeability")->SetMesh(S->GetMesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList wrm_plist = plist_.sublist("water retention evaluator");
  Teuchos::RCP<FlowRelations::WRMRichardsEvaluator> wrm =
      Teuchos::rcp(new FlowRelations::WRMRichardsEvaluator(wrm_plist));
  S->SetFieldEvaluator("saturation_liquid", wrm);
  S->SetFieldEvaluator("saturation_gas", wrm);

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



// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void Richards::initialize(const Teuchos::Ptr<State>& S) {
  // initialize BDF stuff and physical domain stuff
  PKPhysicalBDFBase::initialize(S);

  // initialize boundary conditions
  int nfaces = S->GetMesh()->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::MFD_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  // Initialize face values as the mean of neighboring cell values.
  Teuchos::RCP<CompositeVector> pres = S->GetFieldData(key_, name_);
  DeriveFaceValuesFromCellValues_(pres);

  // Set extra fields as initialized -- these don't currently have evaluators.
  S->GetFieldData("numerical_rel_perm",name_)->PutScalar(1.0);
  S->GetField("numerical_rel_perm",name_)->set_initialized();

  S->GetFieldData("darcy_flux", name_)->PutScalar(0.0);
  S->GetField("darcy_flux", name_)->set_initialized();

  S->GetFieldData("darcy_velocity", name_)->PutScalar(0.0);
  S->GetField("darcy_velocity", name_)->set_initialized();

  // absolute perm
  SetAbsolutePermeabilityTensor_(S);

  // operators
  matrix_->CreateMFDmassMatrices(K_.ptr());
  preconditioner_->CreateMFDmassMatrices(K_.ptr());
};


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void Richards::commit_state(double dt, const Teuchos::RCP<State>& S) {
  // update the rel perm on cells.
  bool update = S->GetFieldEvaluator("relative_permeability")->HasFieldChanged(S.ptr(), name_);
  if (update) UpdatePermeabilityData_(S);
  Teuchos::RCP<const CompositeVector> rel_perm = S->GetFieldData("numerical_rel_perm", name_);

  // update the flux
  Teuchos::RCP<const CompositeVector> pres =
    S->GetFieldData(key_);
  Teuchos::RCP<CompositeVector> darcy_flux =
    S->GetFieldData("darcy_flux", name_);

  matrix_->CreateMFDstiffnessMatrices(*rel_perm);
  matrix_->DeriveFlux(*pres, darcy_flux);
  AddGravityFluxesToVector_(S, darcy_flux);
};


// -----------------------------------------------------------------------------
// Update any diagnostic variables prior to vis (in this case velocity field).
// -----------------------------------------------------------------------------
void Richards::calculate_diagnostics(const Teuchos::RCP<State>& S) {
  // update the cell velocities
  Teuchos::RCP<CompositeVector> velocity = S->GetFieldData("darcy_velocity", name_);
  Teuchos::RCP<const CompositeVector> flux = S->GetFieldData("darcy_flux");
  matrix_->DeriveCellVelocity(*flux, velocity);
};


// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
void Richards::UpdatePermeabilityData_(const Teuchos::RCP<State>& S) {
  Teuchos::RCP<const CompositeVector> rel_perm = S->GetFieldData("relative_permeability");
  //  std::cout << "REL PERM:" << std::endl;
  //  rel_perm->Print(std::cout);

  // get the upwinding done
  upwinding_->Update(S.ptr());

  // patch up the BCs
  //  Teuchos::RCP<const CompositeVector> rel_perm = S->GetFieldData("relative_permeability");
  Teuchos::RCP<CompositeVector> num_rel_perm = S->GetFieldData("numerical_rel_perm", name_);


  for (int c=0; c!=num_rel_perm->size("cell",false); ++c) {
    if (boost::math::isnan<double>((*num_rel_perm)("cell",c))) {
      std::cout << "NaN in cell rel perm." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }
  }
  for (int f=0; f!=num_rel_perm->size("face",false); ++f) {
    if (boost::math::isnan<double>((*num_rel_perm)("face",f))) {
      std::cout << "NaN in face rel perm." << std::endl;
      Errors::Message m("Cut time step");
      Exceptions::amanzi_throw(m);
    }
  }


  for (int f=0; f!=num_rel_perm->size("face"); ++f) {
    AmanziMesh::Entity_ID_List cells;
    num_rel_perm->mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
    if (cells.size() < 2) {
      // just grab the cell inside's perm... this will need to be fixed eventually.
      (*num_rel_perm)("face",f) = (*rel_perm)("cell",cells[0]);
    }
  }

  // scale cells by n/visc
  S->GetFieldEvaluator("molar_density_liquid")->HasFieldChanged(S.ptr(), name_);
  S->GetFieldEvaluator("viscosity_liquid")->HasFieldChanged(S.ptr(), name_);
  Teuchos::RCP<const CompositeVector> n_liq = S->GetFieldData("molar_density_liquid");
  Teuchos::RCP<const CompositeVector> visc = S->GetFieldData("viscosity_liquid");
  int ncells = num_rel_perm->size("cell");
  for (int c=0; c!=ncells; ++c) {
    (*num_rel_perm)("cell",c) *= (*n_liq)("cell",c) / (*visc)("cell",c);
  }

  //  num_rel_perm->Print(std::cout);
};


// -----------------------------------------------------------------------------
// Interpolate pressure ICs on cells to ICs for lambda (faces).
// -----------------------------------------------------------------------------
void Richards::DeriveFaceValuesFromCellValues_(const Teuchos::RCP<CompositeVector>& pres) {
  AmanziMesh::Entity_ID_List cells;
  pres->ScatterMasterToGhosted("cell");

  int f_owned = pres->size("face");
  for (int f=0; f!=f_owned; ++f) {
    cells.clear();
    pres->mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    double face_value = 0.0;
    for (int n=0; n!=ncells; ++n) {
      face_value += (*pres)("cell",cells[n]);
    }
    (*pres)("face",f) = face_value / ncells;
  }
};


// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void Richards::UpdateBoundaryConditions_() {
  for (int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::MFD_BC_NULL;
    bc_values_[n] = 0.0;
  }

  Functions::BoundaryFunction::Iterator bc;
  for (bc=bc_pressure_->begin(); bc!=bc_pressure_->end(); ++bc) {
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
void
Richards::ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& pres) {
  int nfaces = pres->size("face");
  for (int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] == Operators::MFD_BC_DIRICHLET) {
      (*pres)("face",f) = bc_values_[f];
    }
  }
};


} // namespace
} // namespace
