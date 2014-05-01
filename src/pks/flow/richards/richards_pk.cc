/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code.
License: BSD
Authors: Neil Carlson (version 1)
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Ethan Coon (ATS version) (ecoon@lanl.gov)
------------------------------------------------------------------------- */
#include "boost/math/special_functions/fpclassify.hpp"

#include "Epetra_Import.h"

#include "flow_bc_factory.hh"

#include "Point.hh"

#include "upwind_cell_centered.hh"
#include "upwind_arithmetic_mean.hh"
#include "upwind_total_flux.hh"
#include "upwind_gravity_flux.hh"

#include "composite_vector_function.hh"
#include "composite_vector_function_factory.hh"

#include "MatrixMFD_Factory.hh"

#include "predictor_delegate_bc_flux.hh"
#include "wrm_evaluator.hh"
#include "rel_perm_evaluator.hh"
#include "richards_water_content.hh"

#include "richards.hh"

#define DEBUG_RES_FLAG 0


namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
Richards::Richards(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                   Teuchos::ParameterList& FElist,
                   const Teuchos::RCP<TreeVector>& solution) :
    PKDefaultBase(plist, FElist, solution),
    PKPhysicalBDFBase(plist, FElist, solution),
    coupled_to_surface_via_head_(false),
    coupled_to_surface_via_flux_(false),
    infiltrate_only_if_unfrozen_(false),
    modify_predictor_with_consistent_faces_(false),
    modify_predictor_wc_(false),
    modify_predictor_bc_flux_(false),
    upwind_from_prev_flux_(false),
    precon_wc_(false),
    niter_(0),
    dynamic_mesh_(false),
    clobber_surf_kr_(false),
    vapor_diffusion_(false),
    perm_scale_(1.)
{
  // set a few parameters before setup
  plist_->set("primary variable key", "pressure");
  plist_->sublist("primary variable evaluator").set("manage communication", true);
}

// -------------------------------------------------------------
// Setup data
// -------------------------------------------------------------
void Richards::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBDFBase::setup(S);
  SetupRichardsFlow_(S);
  SetupPhysicalEvaluators_(S);

  flux_tol_ = plist_->get<double>("flux tolerance", 1.e-4);
};


// -------------------------------------------------------------
// Pieces of the construction process that are common to all
// Richards-like PKs.
// -------------------------------------------------------------
void Richards::SetupRichardsFlow_(const Teuchos::Ptr<State>& S) {

  // Require fields and evaluators for those fields.
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  // -- primary variable: pressure on both cells and faces, ghosted, with 1 dof
  S->RequireField(key_, name_)->SetMesh(mesh_)->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);

#if DEBUG_RES_FLAG
  // -- residuals of various iterations for debugging
  for (unsigned int i=1; i!=23; ++i) {
    std::stringstream namestream;
    namestream << "flow_residual_" << i;
    std::stringstream solnstream;
    solnstream << "flow_solution_" << i;
    S->RequireField(namestream.str(), name_)->SetMesh(mesh_)->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);
    S->RequireField(solnstream.str(), name_)->SetMesh(mesh_)->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);
  }
#endif

  // -- secondary variables, no evaluator used
  S->RequireField("darcy_flux_direction", name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("face", AmanziMesh::FACE, 1);
  S->RequireField("darcy_flux", name_)->SetMesh(mesh_)->SetGhosted()
                                ->SetComponent("face", AmanziMesh::FACE, 1);
  S->RequireField("darcy_velocity", name_)->SetMesh(mesh_)->SetGhosted()
                                ->SetComponent("cell", AmanziMesh::CELL, 3);

  // Get data for non-field quanitites.
  S->RequireFieldEvaluator("cell_volume");
  S->RequireGravity();
  S->RequireScalar("atmospheric_pressure");

  // Create the absolute permeability tensor.
  unsigned int c_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  K_ = Teuchos::rcp(new std::vector<WhetStone::Tensor>(c_owned));
  for (unsigned int c=0; c!=c_owned; ++c) {
    (*K_)[c].init(mesh_->space_dimension(),1);
  }
  // scaling for permeability
  perm_scale_ = plist_->get<double>("permeability rescaling", 1.0);

  // source terms
  is_source_term_ = plist_->get<bool>("source term", false);
  if (is_source_term_) {
    explicit_source_ = plist_->get<bool>("explicit source term", false);
    S->RequireField("mass_source")->SetMesh(mesh_)
        ->AddComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator("mass_source");
  }
  
  // Create the boundary condition data structures.
  Teuchos::ParameterList bc_plist = plist_->sublist("boundary conditions", true);
  FlowBCFactory bc_factory(mesh_, bc_plist);
  bc_pressure_ = bc_factory.CreatePressure();
  bc_flux_ = bc_factory.CreateMassFlux();
  infiltrate_only_if_unfrozen_ = bc_plist.get<bool>("infiltrate only if unfrozen",false);
  bc_seepage_ = bc_factory.CreateSeepageFace();
  bc_seepage_->Compute(0.); // compute at t=0 to set up

  // how often to update the fluxes?
  std::string updatestring = plist_->get<std::string>("update flux mode", "iteration");
  if (updatestring == "iteration") {
    update_flux_ = UPDATE_FLUX_ITERATION;
  } else if (updatestring == "timestep") {
    update_flux_ = UPDATE_FLUX_TIMESTEP;
  } else if (updatestring == "vis") {
    update_flux_ = UPDATE_FLUX_VIS;
  } else if (updatestring == "never") {
    update_flux_ = UPDATE_FLUX_NEVER;
  } else {
    Errors::Message message(std::string("Unknown frequence for updating the overland flux: ")+updatestring);
    Exceptions::amanzi_throw(message);
  }

  // coupling
  // -- coupling done by a Neumann condition
  coupled_to_surface_via_flux_ = plist_->get<bool>("coupled to surface via flux", false);
  if (coupled_to_surface_via_flux_) {
    S->RequireField("surface_subsurface_flux")
        ->SetMesh(S->GetMesh("surface"))
        ->AddComponent("cell", AmanziMesh::CELL, 1);
  }

  // -- coupling done by a Dirichlet condition
  coupled_to_surface_via_head_ = plist_->get<bool>("coupled to surface via head", false);
  if (coupled_to_surface_via_head_) {
    S->RequireField("surface_pressure");
    // override the flux update -- must happen every iteration
    update_flux_ = UPDATE_FLUX_ITERATION;
  }

  // -- Make sure coupling isn't flagged multiple ways.
  ASSERT(!(coupled_to_surface_via_flux_ && coupled_to_surface_via_head_));


  // Create the upwinding method
  S->RequireField("sc_numerical_rel_perm", name_)->SetMesh(mesh_)->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);
  S->GetField("sc_numerical_rel_perm",name_)->set_io_vis(false);



  S->RequireField("numerical_rel_perm", name_)->SetMesh(mesh_)->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);
  S->GetField("numerical_rel_perm",name_)->set_io_vis(false);

  clobber_surf_kr_ = plist_->get<bool>("clobber surface rel perm", false);
  std::string method_name = plist_->get<std::string>("relative permeability method", "upwind with gravity");
  symmetric_ = false;
  if (method_name == "upwind with gravity") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindGravityFlux(name_,
            "relative_permeability", "numerical_rel_perm", K_));
    Krel_method_ = Operators::UPWIND_METHOD_GRAVITY;
  } else if (method_name == "cell centered") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindCellCentered(name_,
            "relative_permeability", "numerical_rel_perm"));
    symmetric_ = true;
    Krel_method_ = Operators::UPWIND_METHOD_CENTERED;
  } else if (method_name == "upwind with Darcy flux") {
    upwind_from_prev_flux_ = plist_->get<bool>("upwind flux from previous iteration", false);
    if (upwind_from_prev_flux_) {
      upwinding_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
                    "relative_permeability", "numerical_rel_perm", "darcy_flux", 1.e-8));
    } else {
      upwinding_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
                    "relative_permeability", "numerical_rel_perm", "darcy_flux_direction", 1.e-8));
    }
    Krel_method_ = Operators::UPWIND_METHOD_TOTAL_FLUX;
  } else if (method_name == "arithmetic mean") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindArithmeticMean(name_,
            "relative_permeability", "numerical_rel_perm"));
    Krel_method_ = Operators::UPWIND_METHOD_ARITHMETIC_MEAN;
  } else {
    std::stringstream messagestream;
    messagestream << "Richards FLow PK has no upwinding method named: " << method_name;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }


  vapor_diffusion_ = plist_->get<bool>("include vapor diffusion", false);
  if (vapor_diffusion_){
    // Create the vapor diffusion vectors
    S->RequireField("vapor_diffusion_pressure", name_)->SetMesh(mesh_)->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
    S->GetField("vapor_diffusion_pressure",name_)->set_io_vis(true);


    S->RequireField("vapor_diffusion_temperature", name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->GetField("vapor_diffusion_temperature",name_)->set_io_vis(true);
  }

  // operator for the diffusion terms
  Teuchos::ParameterList mfd_plist = plist_->sublist("Diffusion");
  scaled_constraint_ = mfd_plist.get<bool>("scaled constraint equation", false);
  matrix_ = Operators::CreateMatrixMFD(mfd_plist, mesh_);
  symmetric_ = false;
  matrix_->set_symmetric(symmetric_);
  matrix_->SymbolicAssembleGlobalMatrices();
  matrix_->InitPreconditioner();

  if (vapor_diffusion_){
    // operator for the vapor diffusion terms
    matrix_vapor_ = Operators::CreateMatrixMFD(mfd_plist, mesh_);
    symmetric_ = false;
    matrix_vapor_ ->set_symmetric(symmetric_);
    matrix_vapor_ ->SymbolicAssembleGlobalMatrices();
    matrix_vapor_ ->InitPreconditioner();
  }

  // operator with no krel for flux direction, consistent faces
  face_matrix_ = Operators::CreateMatrixMFD(mfd_plist, mesh_);
  face_matrix_->set_symmetric(symmetric_);
  face_matrix_->SymbolicAssembleGlobalMatrices();
  face_matrix_->InitPreconditioner();

  // preconditioner for the NKA system
  Teuchos::ParameterList mfd_pc_plist = plist_->sublist("Diffusion PC");
  if (scaled_constraint_ && !mfd_pc_plist.isParameter("scaled constraint equation"))
    mfd_pc_plist.set("scaled constraint equation", scaled_constraint_);
  mfd_preconditioner_ = Operators::CreateMatrixMFD(mfd_pc_plist, mesh_);
  mfd_preconditioner_->set_symmetric(symmetric_);
  mfd_preconditioner_->SymbolicAssembleGlobalMatrices();
  precon_used_ = mfd_pc_plist.isSublist("preconditioner");
  if (precon_used_)  mfd_preconditioner_->InitPreconditioner();

  // wc preconditioner
  precon_wc_ = plist_->get<bool>("precondition using WC", false);

  // predictors for time integration
  modify_predictor_with_consistent_faces_ =
    plist_->get<bool>("modify predictor with consistent faces", false);
  modify_predictor_bc_flux_ =
    plist_->get<bool>("modify predictor for flux BCs", false);
  modify_predictor_first_bc_flux_ =
    plist_->get<bool>("modify predictor for initial flux BCs", false);
  modify_predictor_wc_ =
    plist_->get<bool>("modify predictor via water content", false);
}


// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void Richards::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  // -- Absolute permeability.
  //       For now, we assume scalar permeability.  This will change.
  S->RequireField("permeability")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("permeability");

  // -- water content, and evaluator
  S->RequireField("water_content")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("water_content");

  // -- Water retention evaluators
  // -- saturation
  Teuchos::ParameterList wrm_plist = plist_->sublist("water retention evaluator");
  Teuchos::RCP<FlowRelations::WRMEvaluator> wrm =
      Teuchos::rcp(new FlowRelations::WRMEvaluator(wrm_plist));
  S->SetFieldEvaluator("saturation_liquid", wrm);
  S->SetFieldEvaluator("saturation_gas", wrm);

  // -- rel perm
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::BOUNDARY_FACE;
  names2[0] = "cell";
  names2[1] = "boundary_face";

  S->RequireField("relative_permeability")->SetMesh(mesh_)->SetGhosted()
      ->AddComponents(names2, locations2, num_dofs2);
  wrm_plist.set<double>("permeability rescaling", perm_scale_);
  Teuchos::RCP<FlowRelations::RelPermEvaluator> rel_perm_evaluator =
      Teuchos::rcp(new FlowRelations::RelPermEvaluator(wrm_plist, wrm->get_WRMs()));
  S->SetFieldEvaluator("relative_permeability", rel_perm_evaluator);
  wrms_ = wrm->get_WRMs();

  // -- Liquid density and viscosity for the transmissivity.
  S->RequireField("molar_density_liquid")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("molar_density_liquid");

  // -- liquid mass density for the gravity fluxes
  S->RequireField("mass_density_liquid")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("mass_density_liquid"); // simply picks up the molar density one.
}


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void Richards::initialize(const Teuchos::Ptr<State>& S) {
  // Initialize BDF stuff and physical domain stuff.
  PKPhysicalBDFBase::initialize(S);

  // debugggin cruft
#if DEBUG_RES_FLAG
  for (unsigned int i=1; i!=23; ++i) {
    std::stringstream namestream;
    namestream << "flow_residual_" << i;
    S->GetFieldData(namestream.str(),name_)->PutScalar(0.);
    S->GetField(namestream.str(),name_)->set_initialized();

    std::stringstream solnstream;
    solnstream << "flow_solution_" << i;
    S->GetFieldData(solnstream.str(),name_)->PutScalar(0.);
    S->GetField(solnstream.str(),name_)->set_initialized();
  }
#endif

  // check whether this is a dynamic mesh problem
  if (S->HasField("vertex coordinate")) dynamic_mesh_ = true;

  // Initialize boundary conditions.
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::MATRIX_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  // Set extra fields as initialized -- these don't currently have evaluators,
  // and will be initialized in the call to commit_state()
  S->GetFieldData("numerical_rel_perm",name_)->PutScalar(1.0);
  S->GetField("numerical_rel_perm",name_)->set_initialized();
  S->GetFieldData("sc_numerical_rel_perm",name_)->PutScalar(1.0);
  S->GetField("sc_numerical_rel_perm",name_)->set_initialized();

  if (vapor_diffusion_){
    S->GetFieldData("vapor_diffusion_pressure",name_)->PutScalar(1.0);
    S->GetField("vapor_diffusion_pressure",name_)->set_initialized();
    S->GetFieldData("vapor_diffusion_temperature",name_)->PutScalar(1.0);
    S->GetField("vapor_diffusion_temperature",name_)->set_initialized();
  }


  S->GetFieldData("darcy_flux", name_)->PutScalar(0.0);
  S->GetField("darcy_flux", name_)->set_initialized();
  S->GetFieldData("darcy_flux_direction", name_)->PutScalar(0.0);
  S->GetField("darcy_flux_direction", name_)->set_initialized();
  S->GetFieldData("darcy_velocity", name_)->PutScalar(0.0);
  S->GetField("darcy_velocity", name_)->set_initialized();

  // absolute perm
  SetAbsolutePermeabilityTensor_(S);

  // operators
  matrix_->CreateMFDmassMatrices(K_.ptr());
  mfd_preconditioner_->CreateMFDmassMatrices(K_.ptr());

  if (vapor_diffusion_){
    //vapor diffusion
    matrix_vapor_->CreateMFDmassMatrices(Teuchos::null);
    // residual vector for vapor diffusion
    res_vapor = Teuchos::rcp(new CompositeVector(*S->GetFieldData("pressure"))); 
  }
  

  face_matrix_->CreateMFDmassMatrices(K_.ptr());
  face_matrix_->CreateMFDstiffnessMatrices(Teuchos::null);


};


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void Richards::commit_state(double dt, const Teuchos::RCP<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Commiting state." << std::endl;

  niter_ = 0;

  bool update = UpdatePermeabilityData_(S.ptr());
  update |= S->GetFieldEvaluator(key_)->HasFieldChanged(S.ptr(), name_);
  update |= S->GetFieldEvaluator("mass_density_liquid")->HasFieldChanged(S.ptr(), name_);

  if (update_flux_ == UPDATE_FLUX_TIMESTEP ||
      (update_flux_ == UPDATE_FLUX_ITERATION && update)) {

    Teuchos::RCP<const CompositeVector> rel_perm =
      S->GetFieldData("numerical_rel_perm");
    // update the stiffness matrix
    matrix_->CreateMFDstiffnessMatrices(rel_perm.ptr());

    // derive fluxes
    Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
    Teuchos::RCP<const CompositeVector> rho = S->GetFieldData("mass_density_liquid");
    Teuchos::RCP<const Epetra_Vector> gvec = S->GetConstantVectorData("gravity");
    Teuchos::RCP<CompositeVector> flux = S->GetFieldData("darcy_flux", name_);
    matrix_->DeriveFlux(*pres, flux.ptr());
    AddGravityFluxesToVector_(gvec.ptr(), rel_perm.ptr(), rho.ptr(), flux.ptr());
  }

  // As a diagnostic, calculate the mass balance error
#if DEBUG_FLAG
  if (S_next_ != Teuchos::null) {
    Teuchos::RCP<const CompositeVector> wc1 = S_next_->GetFieldData("water_content");
    Teuchos::RCP<const CompositeVector> wc0 = S_->GetFieldData("water_content");
    Teuchos::RCP<const CompositeVector> darcy_flux = S->GetFieldData("darcy_flux", name_);
    CompositeVector error(*wc1);

    for (unsigned int c=0; c!=error.size("cell"); ++c) {
      error("cell",c) = (*wc1)("cell",c) - (*wc0)("cell",c);

      AmanziMesh::Entity_ID_List faces;
      std::vector<int> dirs;
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      for (unsigned int n=0; n!=faces.size(); ++n) {
        error("cell",c) += (*darcy_flux)("face",faces[n]) * dirs[n] * dt;
      }
    }

    double einf(0.0);
    error.NormInf(&einf);

    // VerboseObject stuff.
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "Final Mass Balance Error: " << einf << std::endl;
  }
#endif
};


// -----------------------------------------------------------------------------
// Update any diagnostic variables prior to vis (in this case velocity field).
// -----------------------------------------------------------------------------
void Richards::calculate_diagnostics(const Teuchos::RCP<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Calculating diagnostic variables." << std::endl;

  // update the cell velocities
  if (update_flux_ == UPDATE_FLUX_VIS) {
    Teuchos::RCP<const CompositeVector> rel_perm =
      S->GetFieldData("numerical_rel_perm");
    // update the stiffness matrix
    matrix_->CreateMFDstiffnessMatrices(rel_perm.ptr());

    // derive fluxes
    Teuchos::RCP<CompositeVector> flux = S->GetFieldData("darcy_flux", name_);
    Teuchos::RCP<const CompositeVector> pres = S->GetFieldData("pressure");
    Teuchos::RCP<const CompositeVector> rho = S->GetFieldData("mass_density_liquid");
    Teuchos::RCP<const Epetra_Vector> gvec = S->GetConstantVectorData("gravity");
    matrix_->DeriveFlux(*pres, flux.ptr());
    AddGravityFluxesToVector_(gvec.ptr(), rel_perm.ptr(), rho.ptr(), flux.ptr());
  }

  if (update_flux_ != UPDATE_FLUX_NEVER) {
    Teuchos::RCP<CompositeVector> darcy_velocity = S->GetFieldData("darcy_velocity", name_);
    Teuchos::RCP<const CompositeVector> flux = S->GetFieldData("darcy_flux");
    matrix_->DeriveCellVelocity(*flux, darcy_velocity.ptr());

    S->GetFieldEvaluator("molar_density_liquid")->HasFieldChanged(S.ptr(), name_);
    const Epetra_MultiVector& nliq_c = *S->GetFieldData("molar_density_liquid")
        ->ViewComponent("cell",false);

    Epetra_MultiVector& vel_c = *darcy_velocity->ViewComponent("cell",false);
    unsigned int ncells = vel_c.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      for (int n=0; n!=vel_c.NumVectors(); ++n) {
        vel_c[n][c] /= nliq_c[0][c];
      }
    }

  }
};


// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
bool Richards::UpdatePermeabilityData_(const Teuchos::Ptr<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating permeability?";

  Teuchos::RCP<CompositeVector> sc_uw_rel_perm = S->GetFieldData("sc_numerical_rel_perm", name_);
  Teuchos::RCP<CompositeVector> uw_rel_perm = S->GetFieldData("numerical_rel_perm", name_);
  Teuchos::RCP<const CompositeVector> rel_perm = S->GetFieldData("relative_permeability");
  bool update_perm = S->GetFieldEvaluator("relative_permeability")
      ->HasFieldChanged(S, name_);

  // requirements due to the upwinding method
  if (Krel_method_ == Operators::UPWIND_METHOD_TOTAL_FLUX) {
    bool update_dir = S->GetFieldEvaluator("mass_density_liquid")
        ->HasFieldChanged(S, name_);
    update_dir |= S->GetFieldEvaluator(key_)->HasFieldChanged(S, name_);

    if (update_dir) {
      // update the direction of the flux -- note this is NOT the flux
      Teuchos::RCP<CompositeVector> flux_dir =
          S->GetFieldData("darcy_flux_direction", name_);

      // Derive the pressure fluxes
      Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(key_);
      face_matrix_->DeriveFlux(*pres, flux_dir.ptr());

      // Add in the gravity fluxes
      Teuchos::RCP<const Epetra_Vector> gvec = S->GetConstantVectorData("gravity");
      Teuchos::RCP<const CompositeVector> rho = S->GetFieldData("mass_density_liquid");
      // AddGravityFluxesToVector_(gvec.ptr(), sc_uw_rel_perm.ptr(), rho.ptr(), flux_dir.ptr());
      AddGravityFluxesToVector_(gvec.ptr(), Teuchos::null, rho.ptr(), flux_dir.ptr());
    }

    update_perm |= update_dir;
  }

  if (update_perm) {
    // Move rel perm on boundary_faces into uw_rel_perm on faces
    const Epetra_Import& vandelay = mesh_->exterior_face_importer();
    const Epetra_MultiVector& rel_perm_bf =
        *rel_perm->ViewComponent("boundary_face",false);
    {
      Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->ViewComponent("face",false);
      uw_rel_perm_f.Export(rel_perm_bf, vandelay, Insert);
    }

    // Upwind, only overwriting boundary faces if the wind says to do so.
    upwinding_->Update(S);

    if (clobber_surf_kr_) {
      Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->ViewComponent("face",false);
      uw_rel_perm_f.Export(rel_perm_bf, vandelay, Insert);
    }
  }

  // debugging
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << " " << update_perm << std::endl;
  }
  return update_perm;
};


// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void Richards::UpdateBoundaryConditions_() {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating BCs." << std::endl;

  for (unsigned int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::MATRIX_BC_NULL;
    bc_values_[n] = 0.0;
  }

  // Dirichlet boundary conditions
  Functions::BoundaryFunction::Iterator bc;
  for (bc=bc_pressure_->begin(); bc!=bc_pressure_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
    bc_values_[f] = bc->second;
  }

  if (!infiltrate_only_if_unfrozen_) {
    // Standard Neuman boundary conditions
    for (bc=bc_flux_->begin(); bc!=bc_flux_->end(); ++bc) {
      int f = bc->first;
      bc_markers_[f] = Operators::MATRIX_BC_FLUX;
      bc_values_[f] = bc->second;
    }
  } else {
    // Neumann boundary conditions that turn off if temp < freezing
    const Epetra_MultiVector& temp = *S_next_->GetFieldData("temperature")->ViewComponent("face");
    for (bc=bc_flux_->begin(); bc!=bc_flux_->end(); ++bc) {
      int f = bc->first;
      bc_markers_[f] = Operators::MATRIX_BC_FLUX;
      if (temp[0][f] > 273.15) {
        bc_values_[f] = bc->second;
      } else {
        bc_values_[f] = 0.;
      }
    }
  }

  // seepage face -- pressure <= p_atm, outward mass flux >= 0
  const Epetra_MultiVector pressure = *S_next_->GetFieldData(key_)->ViewComponent("face");
  const double& p_atm = *S_next_->GetScalarData("atmospheric_pressure");
  for (bc=bc_seepage_->begin(); bc!=bc_seepage_->end(); ++bc) {
    int f = bc->first;
    if (pressure[0][f] < p_atm) {
      bc_markers_[f] = Operators::MATRIX_BC_FLUX;
      bc_values_[f] = bc->second;
    } else {
      bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
      bc_values_[f] = p_atm;
    }
  }

  // surface coupling
  if (coupled_to_surface_via_head_) {
    // Face is Dirichlet with value of surface head
    Teuchos::RCP<const AmanziMesh::Mesh> surface = S_next_->GetMesh("surface");
    const Epetra_MultiVector& head = *S_next_->GetFieldData("surface_pressure")
        ->ViewComponent("cell",false);

    unsigned int ncells_surface = head.MyLength();
    for (unsigned int c=0; c!=ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face
      AmanziMesh::Entity_ID f =
        surface->entity_get_parent(AmanziMesh::CELL, c);

      // -- set that value to dirichlet
      bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
      bc_values_[f] = head[0][c];
    }
  }

  // surface coupling
  if (coupled_to_surface_via_flux_) {
    // Face is Neumann with value of surface residual
    Teuchos::RCP<const AmanziMesh::Mesh> surface = S_next_->GetMesh("surface");
    const Epetra_MultiVector& flux = *S_next_->GetFieldData("surface_subsurface_flux")
        ->ViewComponent("cell",false);

    unsigned int ncells_surface = flux.MyLength();
    for (unsigned int c=0; c!=ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face
      AmanziMesh::Entity_ID f =
        surface->entity_get_parent(AmanziMesh::CELL, c);

      // -- set that value to Neumann
      bc_markers_[f] = Operators::MATRIX_BC_FLUX;
      bc_values_[f] = flux[0][c] / mesh_->face_area(f);
      // NOTE: flux[0][c] is in units of mols / s, where as Neumann BCs are in
      //       units of mols / s / A.  The right A must be chosen, as it is
      //       the subsurface mesh's face area, not the surface mesh's cell
      //       area.
    }
  }
};


// -----------------------------------------------------------------------------
// Add a boundary marker to owned faces.
// -----------------------------------------------------------------------------
void
Richards::ApplyBoundaryConditions_(const Teuchos::Ptr<CompositeVector>& pres) {
  Epetra_MultiVector& pres_f = *pres->ViewComponent("face",false);
  unsigned int nfaces = pres_f.MyLength();
  for (unsigned int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] == Operators::MATRIX_BC_DIRICHLET) {
      pres_f[0][f] = bc_values_[f];
    }
  }
};


bool Richards::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
        Teuchos::RCP<TreeVector> u) {
  
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Modifying predictor:" << std::endl;

  bool changed(false);
  if (modify_predictor_bc_flux_ ||
      (modify_predictor_first_bc_flux_ && S_next_->cycle() == 0)) {
    changed |= ModifyPredictorFluxBCs_(h,u);
  }

  if (modify_predictor_wc_) {
    changed |= ModifyPredictorWC_(h,u);
  }

  if (modify_predictor_with_consistent_faces_) {
    changed |= ModifyPredictorConsistentFaces_(h,u);
  }
  return changed;
}

bool Richards::ModifyPredictorFluxBCs_(double h, Teuchos::RCP<TreeVector> u) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  modifications to deal with nonlinearity at flux BCs" << std::endl;

  if (flux_predictor_ == Teuchos::null) {
    flux_predictor_ = Teuchos::rcp(new PredictorDelegateBCFlux(S_next_, mesh_, matrix_,
            wrms_, &bc_markers_, &bc_values_));
  }

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  UpdatePermeabilityData_(S_next_.ptr());
  Teuchos::RCP<const CompositeVector> rel_perm =
    S_next_->GetFieldData("numerical_rel_perm");
  matrix_->CreateMFDstiffnessMatrices(rel_perm.ptr());
  matrix_->CreateMFDrhsVectors();
  Teuchos::RCP<const CompositeVector> rho = S_next_->GetFieldData("mass_density_liquid");
  Teuchos::RCP<const Epetra_Vector> gvec = S_next_->GetConstantVectorData("gravity");
  AddGravityFluxes_(gvec.ptr(), rel_perm.ptr(), rho.ptr(), matrix_.ptr());
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);

  flux_predictor_->ModifyPredictor(h, u);
  ChangedSolution(); // mark the solution as changed, as modifying with
                      // consistent faces will then get the updated boundary
                      // conditions
  return true;
}

bool Richards::ModifyPredictorConsistentFaces_(double h, Teuchos::RCP<TreeVector> u) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  modifications for consistent face pressures." << std::endl;
  CalculateConsistentFaces(u->Data().ptr());
  return true;
}

bool Richards::ModifyPredictorWC_(double h, Teuchos::RCP<TreeVector> u) {
  ASSERT(0);
  return false;
}


void Richards::CalculateConsistentFacesForInfiltration_(
    const Teuchos::Ptr<CompositeVector>& u) {
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  modifications to deal with nonlinearity at flux BCs" << std::endl;

  if (flux_predictor_ == Teuchos::null) {
    flux_predictor_ = Teuchos::rcp(new PredictorDelegateBCFlux(S_next_, mesh_, matrix_,
            wrms_, &bc_markers_, &bc_values_));
  }

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  bool update = UpdatePermeabilityData_(S_next_.ptr());
  Teuchos::RCP<const CompositeVector> rel_perm =
      S_next_->GetFieldData("numerical_rel_perm");
  matrix_->CreateMFDstiffnessMatrices(rel_perm.ptr());
  matrix_->CreateMFDrhsVectors();
  Teuchos::RCP<const CompositeVector> rho = S_next_->GetFieldData("mass_density_liquid");
  Teuchos::RCP<const Epetra_Vector> gvec = S_next_->GetConstantVectorData("gravity");
  AddGravityFluxes_(gvec.ptr(), rel_perm.ptr(), rho.ptr(), matrix_.ptr());
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);

  flux_predictor_->ModifyPredictor(u);
}

void Richards::CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();

  // update the rel perm according to the scheme of choice
  ChangedSolution();
  UpdatePermeabilityData_(S_next_.ptr());

  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_();

  Teuchos::RCP<CompositeVector> rel_perm =
      S_next_->GetFieldData("sc_numerical_rel_perm", name_);

  Teuchos::RCP<const CompositeVector> rho =
      S_next_->GetFieldData("mass_density_liquid");
  Teuchos::RCP<const Epetra_Vector> gvec =
      S_next_->GetConstantVectorData("gravity");

  // Update the preconditioner with darcy and gravity fluxes
  matrix_->CreateMFDstiffnessMatrices(rel_perm.ptr());
  matrix_->CreateMFDrhsVectors();
  AddGravityFluxes_(gvec.ptr(), rel_perm.ptr(), rho.ptr(), matrix_.ptr());

  // skip accumulation terms, they're not needed

  // Assemble
  matrix_->ApplyBoundaryConditions(bc_markers_, bc_values_);
  matrix_->AssembleGlobalMatrices();

  // derive the consistent faces, involves a solve
  matrix_->UpdateConsistentFaceConstraints(u.ptr());
}

// -----------------------------------------------------------------------------
// Check admissibility of the solution guess.
// -----------------------------------------------------------------------------
bool Richards::IsAdmissible(Teuchos::RCP<const TreeVector> up) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Checking admissibility..." << std::endl;

  // For some reason, wandering PKs break most frequently with an unreasonable
  // presure.  This simply tries to catch that before it happens.
  Teuchos::RCP<const CompositeVector> pres = up->Data();

  const Epetra_MultiVector& pres_v = *pres->ViewComponent("cell",false);
  double minp(0.), maxp(0.);
  int ierr = pres_v.MinValue(&minp);
  ierr |= pres_v.MaxValue(&maxp);

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "    Admissible p? (min/max): " << minp << ",  " << maxp << std::endl;
  }

  if (ierr || minp < -1.e8 || maxp > 1.e8) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << " is not admissible, as it is not within bounds of constitutive models: min(p) = " << minp << ", max(p) = " << maxp << std::endl;
    }
    return false;
  }
  return true;
}

} // namespace
} // namespace
