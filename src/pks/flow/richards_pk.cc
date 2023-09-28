/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code.
License: BSD
Authors: Neil Carlson (version 1)
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Ethan Coon (ATS version) (ecoon@lanl.gov)
------------------------------------------------------------------------- */
#include "boost/math/special_functions/fpclassify.hpp"

#include "boost/algorithm/string/predicate.hpp"

#include "Epetra_Import.h"

#include "flow_bc_factory.hh"

#include "Point.hh"

#include "upwind_cell_centered.hh"
#include "upwind_arithmetic_mean.hh"
#include "upwind_total_flux.hh"
#include "upwind_gravity_flux.hh"

#include "CompositeVectorFunction.hh"
#include "CompositeVectorFunctionFactory.hh"

#include "predictor_delegate_bc_flux.hh"
#include "wrm_evaluator.hh"
#include "rel_perm_evaluator.hh"
#include "richards_water_content_evaluator.hh"
#include "OperatorDefs.hh"
#include "BoundaryFlux.hh"
#include "pk_helpers.hh"

#include "richards.hh"

#define DEBUG_RES_FLAG 0


namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------

Richards::Richards(Teuchos::ParameterList& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& glist,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& solution) :
    PK(pk_tree, glist,  S, solution),
    PK_PhysicalBDF_Default(pk_tree, glist,  S, solution),
    coupled_to_surface_via_head_(false),
    coupled_to_surface_via_flux_(false),
    infiltrate_only_if_unfrozen_(false),
    modify_predictor_with_consistent_faces_(false),
    modify_predictor_wc_(false),
    modify_predictor_bc_flux_(false),
    modify_predictor_first_bc_flux_(false),
    upwind_from_prev_flux_(false),
    dynamic_mesh_(false),
    clobber_boundary_flux_dir_(false),
    vapor_diffusion_(false),
    perm_scale_(1.),
    jacobian_(false),
    jacobian_lag_(0),
    iter_(0),
    iter_counter_time_(0.)
{
  // set a default absolute tolerance
  if (!plist_->isParameter("absolute error tolerance"))
    plist_->set("absolute error tolerance", .5 * .1 * 55000.); // phi * s * nl

  // get field names
  conserved_key_ = Keys::readKey(*plist_, domain_, "conserved", "water_content");
  mass_dens_key_ = Keys::readKey(*plist_, domain_, "mass density", "mass_density_liquid");
  molar_dens_key_ = Keys::readKey(*plist_, domain_, "molar density", "molar_density_liquid");
  perm_key_ = Keys::readKey(*plist_, domain_, "permeability", "permeability");
  coef_key_ = Keys::readKey(*plist_, domain_, "conductivity", "relative_permeability");
  uw_coef_key_ = Keys::readKey(*plist_, domain_, "upwinded conductivity", "upwind_relative_permeability");
  flux_key_ = Keys::readKey(*plist_, domain_, "darcy flux", "mass_flux");
  flux_dir_key_ = Keys::readKey(*plist_, domain_, "darcy flux direction", "mass_flux_direction");
  velocity_key_ = Keys::readKey(*plist_, domain_, "darcy velocity", "darcy_velocity");
  sat_key_ = Keys::readKey(*plist_, domain_, "saturation", "saturation_liquid");
  sat_gas_key_ = Keys::readKey(*plist_, domain_, "saturation gas", "saturation_gas");
  sat_ice_key_ = Keys::readKey(*plist_, domain_, "saturation ice", "saturation_ice");

  std::cout << "Richards, domain_ = " << domain_ << std::endl;
  std::cout << "Richards, conserved_key_ = " << conserved_key_ << std::endl;

  // set up an additional primary variable evaluator for flux
  Teuchos::ParameterList& pv_sublist = S->GetEvaluatorList(flux_key_);
  pv_sublist.set("field evaluator type", "primary variable");
}

// -------------------------------------------------------------
// Setup data
// -------------------------------------------------------------
void Richards::Setup(const Teuchos::Ptr<State>& S)
{

  std::cout << "setup richards START" << std::endl;

  PK_PhysicalBDF_Default::Setup(S);
  std::cout << "after PK_PhysicalBDF_Default::Setup(S)" << std::endl;
  SetupRichardsFlow_(S);
  std::cout << "SetupRichardsFlow_(S)" << std::endl;
  SetupPhysicalEvaluators_(S);
  std::cout << "SetupPhysicalEvaluators_(S)" << std::endl;

  std::cout << "setup richards DONE" << std::endl;

};


// -------------------------------------------------------------
// Pieces of the construction process that are common to all
// Richards-like PKs.
// -------------------------------------------------------------
void Richards::SetupRichardsFlow_(const Teuchos::Ptr<State>& S)
{
  // Get data for special-case entities.
  S->RequireField(cell_vol_key_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(cell_vol_key_);
  S->RequireGravity();
  S->RequireScalar("atmospheric_pressure");

  // Set up Operators
  // -- boundary conditions
  Teuchos::ParameterList bc_plist = plist_->sublist("boundary conditions", true);
  FlowBCFactory bc_factory(mesh_, bc_plist);
  bc_pressure_ = bc_factory.CreatePressure();
  bc_head_ = bc_factory.CreateHead();
  bc_level_ = bc_factory.CreateFixedLevel();
  bc_flux_ = bc_factory.CreateMassFlux();
  bc_seepage_ = bc_factory.CreateSeepageFacePressure();
  bc_seepage_->Compute(0.); // compute at t=0 to set up
  std::tie(bc_seepage_infilt_explicit_, bc_seepage_infilt_) =
    bc_factory.CreateSeepageFacePressureWithInfiltration();
  bc_seepage_infilt_->Compute(0.); // compute at t=0 to set up
  bc_rho_water_ = bc_plist.get<double>("hydrostatic water density [kg m^-3]",1000.);

  // scaling for permeability
  perm_scale_ = plist_->get<double>("permeability rescaling", 1.e7);

  // permeability type - scalar or tensor?
  Teuchos::ParameterList& perm_list = S->GetEvaluatorList(perm_key_);
  std::string perm_type = perm_list.get<std::string>("permeability type", "scalar");
  if (perm_type == "scalar") {
    perm_tensor_rank_ = 1;
    num_perm_vals_ = 1;
  } else if (perm_type == "horizontal and vertical") {
    perm_tensor_rank_ = 2;
    num_perm_vals_ = 2;
  } else if (perm_type == "diagonal tensor") {
    perm_tensor_rank_ = 2;
    num_perm_vals_ = mesh_->space_dimension();
  } else if (perm_type == "full tensor") {
    perm_tensor_rank_ = 2;
    num_perm_vals_ = (mesh_->space_dimension() == 3) ? 6 : 3;
  } else {
    Errors::Message message("`permeability type` must be one of the following: \"scalar\", \"diagonal tensor\", \"full tensor\", or \"horizontal and vertical\".");
    Exceptions::amanzi_throw(message);
  }

  // -- linear tensor coefficients
  unsigned int c_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  K_ = Teuchos::rcp(new std::vector<WhetStone::Tensor>(c_owned));
  for (unsigned int c=0; c!=c_owned; ++c) {
    (*K_)[c].Init(mesh_->space_dimension(),perm_tensor_rank_);
  }

  // -- nonlinear coefficients/upwinding
  if (plist_->isParameter("surface rel perm strategy")) {
    clobber_policy_ = plist_->get<std::string>("surface rel perm strategy");
  } else if (plist_->get<bool>("clobber surface rel perm", false)) {
    clobber_policy_ = "clobber";
  } else if (plist_->get<bool>("max surface rel perm", false)) {
    clobber_policy_ = "max";
  } else if (plist_->get<bool>("unsaturated clobber surface rel perm", false)) {
    clobber_policy_ = "unsaturated";
  } else {
    clobber_policy_ = "none";
  }
  clobber_boundary_flux_dir_ = plist_->get<bool>("clobber boundary flux direction for upwinding", false);

  std::string method_name = plist_->get<std::string>("relative permeability method", "upwind with Darcy flux");
  if (method_name == "upwind with gravity") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindGravityFlux(name_,
            coef_key_, uw_coef_key_, K_));
    Krel_method_ = Operators::UPWIND_METHOD_GRAVITY;
  } else if (method_name == "cell centered") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindCellCentered(name_,
            coef_key_, uw_coef_key_));
    Krel_method_ = Operators::UPWIND_METHOD_CENTERED;
  } else if (method_name == "upwind with Darcy flux") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
            coef_key_, uw_coef_key_, flux_dir_key_, 1.e-5));
    Krel_method_ = Operators::UPWIND_METHOD_TOTAL_FLUX;
  } else if (method_name == "arithmetic mean") {
    upwinding_ = Teuchos::rcp(new Operators::UpwindArithmeticMean(name_,
            coef_key_, uw_coef_key_));
    Krel_method_ = Operators::UPWIND_METHOD_ARITHMETIC_MEAN;
  } else {
    std::stringstream messagestream;
    messagestream << "Richards Flow PK has no upwinding method named: " << method_name;
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  // -- require the data on appropriate locations
  std::string coef_location = upwinding_->CoefficientLocation();
  if (coef_location == "upwind: face") {
    S->RequireField(uw_coef_key_, name_)->SetMesh(mesh_)
        ->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 1);
  } else if (coef_location == "standard: cell") {
    S->RequireField(uw_coef_key_, name_)->SetMesh(mesh_)
        ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
  } else {
    Errors::Message message("Unknown upwind coefficient location in Richards flow.");
    Exceptions::amanzi_throw(message);
  }
  S->GetField(uw_coef_key_,name_)->set_io_vis(false);

  // -- create the forward operator for the diffusion term
  Teuchos::ParameterList& mfd_plist = plist_->sublist("diffusion");
  mfd_plist.set("nonlinear coefficient", coef_location);
  mfd_plist.set("gravity", true);

  Operators::PDE_DiffusionFactory opfactory;
  matrix_diff_ = opfactory.CreateWithGravity(mfd_plist, mesh_, bc_);
  matrix_ = matrix_diff_->global_operator();

  // -- create the operator, data for flux directions
  Teuchos::ParameterList face_diff_list(mfd_plist);
  face_diff_list.set("nonlinear coefficient", "none");
  face_matrix_diff_ = opfactory.CreateWithGravity(face_diff_list, mesh_, bc_);

  S->RequireField(flux_dir_key_, name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("face", AmanziMesh::FACE, 1);

  // -- create the operators for the preconditioner
  //    diffusion
  // NOTE: Can this be a clone of the primary operator? --etc
  Teuchos::ParameterList& mfd_pc_plist = plist_->sublist("diffusion preconditioner");
  mfd_pc_plist.set("nonlinear coefficient", coef_location);
  mfd_pc_plist.set("gravity", true);
  if (!mfd_pc_plist.isParameter("discretization primary"))
    mfd_pc_plist.set("discretization primary", mfd_plist.get<std::string>("discretization primary"));
  if (!mfd_pc_plist.isParameter("discretization secondary") && mfd_plist.isParameter("discretization secondary"))
    mfd_pc_plist.set("discretization secondary", mfd_plist.get<std::string>("discretization secondary"));
  if (!mfd_pc_plist.isParameter("schema") && mfd_plist.isParameter("schema"))
    mfd_pc_plist.set("schema", mfd_plist.get<Teuchos::Array<std::string> >("schema"));
  if (mfd_pc_plist.get<bool>("include Newton correction", false)) {
    if (mfd_pc_plist.get<std::string>("discretization primary") == "fv: default") {
      mfd_pc_plist.set("Newton correction", "true Jacobian");
    } else {
      mfd_pc_plist.set("Newton correction", "approximate Jacobian");
    }
  }

  precon_used_ = plist_->isSublist("preconditioner") ||
    plist_->isSublist("inverse") ||
    plist_->isSublist("linear solver");
  if (precon_used_) {
    mfd_pc_plist.set("inverse", plist_->sublist("inverse"));
    // old style... deprecate me!
    mfd_pc_plist.sublist("inverse").setParameters(plist_->sublist("preconditioner"));
    mfd_pc_plist.sublist("inverse").setParameters(plist_->sublist("linear solver"));
  }

  preconditioner_diff_ = opfactory.CreateWithGravity(mfd_pc_plist, mesh_, bc_);
  preconditioner_ = preconditioner_diff_->global_operator();

  //    If using approximate Jacobian for the preconditioner, we also need
  //    derivative information.  For now this means upwinding the derivative.
  jacobian_ = mfd_pc_plist.get<std::string>("Newton correction", "none") != "none";
  if (jacobian_) {
    jacobian_lag_ = mfd_pc_plist.get<int>("Newton correction lag", 0);

    if (mfd_pc_plist.get<std::string>("discretization primary") != "fv: default"){
      // MFD -- upwind required
      dcoef_key_ = Keys::getDerivKey(coef_key_, key_);
      duw_coef_key_ = Keys::getDerivKey(uw_coef_key_, key_);

      S->RequireField(duw_coef_key_, name_)
        ->SetMesh(mesh_)->SetGhosted()
        ->SetComponent("face", AmanziMesh::FACE, 1);

      upwinding_deriv_ = Teuchos::rcp(new Operators::UpwindTotalFlux(name_,
                                      dcoef_key_, duw_coef_key_, flux_dir_key_, 1.e-8));
    } else {
      // FV -- no upwinding
      dcoef_key_ = Keys::getDerivKey(coef_key_, key_);
      duw_coef_key_ = std::string();
    }
  }

  // -- accumulation terms
  Teuchos::ParameterList& acc_pc_plist = plist_->sublist("accumulation preconditioner");
  acc_pc_plist.set<std::string>("entity kind", "cell");
  preconditioner_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(acc_pc_plist, preconditioner_));

  // // -- vapor diffusion terms
  // vapor_diffusion_ = plist_->get<bool>("include vapor diffusion", false);
  // if (vapor_diffusion_){
  //   AMANZI_ASSERT(0); // untested!

  //   // Create the vapor diffusion vectors
  //   S->RequireField("vapor_diffusion_pressure", name_)->SetMesh(mesh_)->SetGhosted()
  //       ->SetComponent("cell", AmanziMesh::CELL, 1);
  //   S->GetField("vapor_diffusion_pressure",name_)->set_io_vis(true);

  //   S->RequireField("vapor_diffusion_temperature", name_)->SetMesh(mesh_)->SetGhosted()
  //     ->SetComponent("cell", AmanziMesh::CELL, 1);
  //   S->GetField("vapor_diffusion_temperature",name_)->set_io_vis(true);

  //   // operator for the vapor diffusion terms
  //   matrix_vapor_ = Operators::CreateMatrixMFD(mfd_plist, mesh_);
  // }

  // source terms
  is_source_term_ = plist_->get<bool>("source term", false);
  if (is_source_term_) {
    if (source_key_.empty())
      source_key_ = Keys::readKey(*plist_, domain_, "source", "water_source");
    source_term_is_differentiable_ =
        plist_->get<bool>("source term is differentiable", true);
    explicit_source_ = plist_->get<bool>("explicit source term", false);
    S->RequireField(source_key_)->SetMesh(mesh_)
        ->AddComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(source_key_);
  }

  // coupling to the surface
  // -- coupling done by a Neumann condition
  coupled_to_surface_via_flux_ = plist_->get<bool>("coupled to surface via flux", false);
  if (coupled_to_surface_via_flux_) {
    std::string domain_surf;
    if (domain_ == "domain" || domain_ == "") {
      domain_surf = plist_->get<std::string>("surface domain name", "surface");
    } else {
      domain_surf = plist_->get<std::string>("surface domain name", "surface_"+domain_);
    }
    ss_flux_key_ = Keys::readKey(*plist_, domain_surf, "surface-subsurface flux", "surface_subsurface_flux");
    S->RequireField(ss_flux_key_)->SetMesh(S->GetMesh(domain_surf))
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  }

  // -- coupling done by a Dirichlet condition
  coupled_to_surface_via_head_ = plist_->get<bool>("coupled to surface via head", false);
  if (coupled_to_surface_via_head_) {
    std::string domain_surf;
    if (domain_ == "domain" || domain_ == "") {
      domain_surf = plist_->get<std::string>("surface domain name", "surface");
    } else {
      domain_surf = plist_->get<std::string>("surface domain name", "surface_"+domain_);
    }
    ss_primary_key_ = Keys::readKey(*plist_, domain_surf, "pressure", "pressure");
    S->RequireField(ss_primary_key_)
      ->SetMesh(S->GetMesh(domain_surf))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }

  // -- Make sure coupling isn't flagged multiple ways.
  if (coupled_to_surface_via_flux_ && coupled_to_surface_via_head_) {
    Errors::Message message("Richards PK requested both flux and head coupling -- choose one.");
    Exceptions::amanzi_throw(message);
  }

  // Require fields and evaluators
  CompositeVectorSpace matrix_cvs = matrix_->RangeMap();
  compute_boundary_values_ = plist_->get<bool>("compute boundary values", false);
  if (compute_boundary_values_)
    matrix_cvs.AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
  S->RequireField(key_, name_)->Update(matrix_cvs)->SetGhosted();

  // -- flux is managed here as a primary variable
  S->RequireField(flux_key_, name_)->SetMesh(mesh_)->SetGhosted()
                                ->SetComponent("face", AmanziMesh::FACE, 1);
  S->RequireFieldEvaluator(flux_key_);

  // -- also need a velocity, but only for vis/diagnostics
  S->RequireField(velocity_key_, name_)->SetMesh(mesh_)->SetGhosted()
                                ->SetComponent("cell", AmanziMesh::CELL, 3);

  // Globalization and other timestep control flags
  // -- predictors
  modify_predictor_with_consistent_faces_ =
    plist_->get<bool>("modify predictor with consistent faces", false);
  modify_predictor_bc_flux_ =
    plist_->get<bool>("modify predictor for flux BCs", false);
  modify_predictor_first_bc_flux_ =
    plist_->get<bool>("modify predictor for initial flux BCs", false);
  modify_predictor_wc_ =
    plist_->get<bool>("modify predictor via water content", false);

  // -- correctors
  p_limit_ = plist_->get<double>("limit correction to pressure change [Pa]", -1.);
  patm_limit_ = plist_->get<double>("limit correction to pressure change when crossing atmospheric [Pa]", -1.);

  // -- valid step controls
  sat_change_limit_ = plist_->get<double>("max valid change in saturation in a time step [-]", -1.);
  sat_ice_change_limit_ = plist_->get<double>("max valid change in ice saturation in a time step [-]", -1.);
}


// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void Richards::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  // -- Absolute permeability.
  //       For now, we assume scalar permeability.  This will change.
  S->RequireField(perm_key_)->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, num_perm_vals_);
  S->RequireFieldEvaluator(perm_key_);

  // -- water content, and evaluator
  S->RequireField(conserved_key_)->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(conserved_key_);
  std::cout << "Richards water content eval conserved_key_ = " << conserved_key_ << std::endl;

  // -- Water retention evaluators
  // This deals with deprecated location for the WRM list (in the PK).  Move it
  // to state.
  if (plist_->isSublist("water retention evaluator")) {
    auto& wrm_plist = S->GetEvaluatorList(sat_key_);
    wrm_plist.setParameters(plist_->sublist("water retention evaluator"));
    wrm_plist.set("field evaluator type", "WRM");
  }
  if (!S->HasFieldEvaluator(coef_key_) && (S->GetEvaluatorList(coef_key_).numParams() == 0)) {
    Teuchos::ParameterList& kr_plist = S->GetEvaluatorList(coef_key_);
    kr_plist.setParameters(S->GetEvaluatorList(sat_key_));
    kr_plist.set("field evaluator type", "WRM rel perm");
  }

  // -- saturation
  S->RequireField(sat_key_)->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField(sat_gas_key_)->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  auto wrm = S->RequireFieldEvaluator(sat_key_);
  S->RequireFieldEvaluator(sat_gas_key_);

  // -- rel perm
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::BOUNDARY_FACE;
  names2[0] = "cell";
  names2[1] = "boundary_face";
  S->RequireField(coef_key_)->SetMesh(mesh_)->SetGhosted()
      ->AddComponents(names2, locations2, num_dofs2);
  S->GetEvaluatorList(coef_key_).set<double>("permeability rescaling", perm_scale_);
  S->RequireFieldEvaluator(coef_key_);

  // -- get the WRM models
  auto wrm_eval = Teuchos::rcp_dynamic_cast<Flow::WRMEvaluator>(wrm);
  AMANZI_ASSERT(wrm_eval != Teuchos::null);
  wrms_ = wrm_eval->get_WRMs();

  // -- Liquid density and viscosity for the transmissivity.
  S->RequireField(molar_dens_key_)->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(molar_dens_key_);

  // -- liquid mass density for the gravity fluxes
  S->RequireField(mass_dens_key_)->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(mass_dens_key_); // simply picks up the molar density one.

}


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void Richards::Initialize(const Teuchos::Ptr<State>& S)
{
  // stash the flux evaluator first -- needs to be done sometime.
  Teuchos::RCP<FieldEvaluator> flux_eval = S->GetFieldEvaluator(flux_key_);
  flux_pvfe_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(flux_eval);
  AMANZI_ASSERT(flux_pvfe_ != Teuchos::null);
  flux_pvfe_->SetFieldAsChanged(S);

  // Initialize via hydrostatic balance
  if (!S->GetField(key_, name_)->initialized())
    InitializeHydrostatic_(S);

  // Initialize BDF stuff and physical domain stuff.
  PK_PhysicalBDF_Default::Initialize(S);

  // debugging cruft
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

  // Set extra fields as initialized -- these don't currently have evaluators,
  // and will be initialized in the call to commit_state()
  S->GetFieldData(uw_coef_key_,name_)->PutScalar(1.0);
  S->GetField(uw_coef_key_,name_)->set_initialized();

  if (!duw_coef_key_.empty()) {
    S->GetFieldData(duw_coef_key_,name_)->PutScalar(1.0);
    S->GetField(duw_coef_key_,name_)->set_initialized();
  }

  S->GetFieldData(flux_key_, name_)->PutScalar(0.0);
  S->GetField(flux_key_, name_)->set_initialized();
  S->GetFieldData(flux_dir_key_, name_)->PutScalar(0.0);
  S->GetField(flux_dir_key_, name_)->set_initialized();
  S->GetFieldData(velocity_key_, name_)->PutScalar(0.0);
  S->GetField(velocity_key_, name_)->set_initialized();

  // absolute perm
  SetAbsolutePermeabilityTensor_(S);

  // operators
  Teuchos::RCP<const Epetra_Vector> gvec = S->GetConstantVectorData("gravity");
  AmanziGeometry::Point g(3);
  g[0] = (*gvec)[0]; g[1] = (*gvec)[1]; g[2] = (*gvec)[2];

  matrix_diff_->SetGravity(g);
  matrix_diff_->SetBCs(bc_, bc_);
  matrix_diff_->SetTensorCoefficient(K_);

  preconditioner_diff_->SetGravity(g);
  preconditioner_diff_->SetBCs(bc_, bc_);
  preconditioner_diff_->SetTensorCoefficient(K_);

  face_matrix_diff_->SetGravity(g);
  face_matrix_diff_->SetBCs(bc_, bc_);
  face_matrix_diff_->SetTensorCoefficient(K_);
  face_matrix_diff_->SetScalarCoefficient(Teuchos::null, Teuchos::null);

  // if (vapor_diffusion_){
  //   //vapor diffusion
  //   matrix_vapor_->CreateMFDmassMatrices(Teuchos::null);
  //   // residual vector for vapor diffusion
  //   res_vapor = Teuchos::rcp(new CompositeVector(*S->GetFieldData(key_)));
  // }
};


void Richards::InitializeHydrostatic_(const Teuchos::Ptr<State>& S)
{
  // constant head over the surface
  if (plist_->sublist("initial condition").isParameter("hydrostatic head [m]")) {
    double head_wt = plist_->sublist("initial condition").get<double>("hydrostatic head [m]");
    double rho = plist_->sublist("initial condition").get<double>("hydrostatic water density [kg m^-3]");
    int ncols = mesh_->num_columns(false);

    Teuchos::RCP<const Epetra_Vector> gvec = S->GetConstantVectorData("gravity");
    int z_index = mesh_->space_dimension() - 1;
    double g = -(*gvec)[z_index];
    double p_atm = *S->GetScalarData("atmospheric_pressure");

    // set pressure on the column of faces and cells
    Teuchos::RCP<CompositeVector> pres = S->GetFieldData(key_, name_);
    Teuchos::RCP<Epetra_IntVector> flags = Teuchos::null;
    if (pres->HasComponent("face"))
      flags = Teuchos::rcp(new Epetra_IntVector(*pres->Map().Map("face",false)));

    {
      Epetra_MultiVector& pres_c = *pres->ViewComponent("cell", false);
      Teuchos::RCP<Epetra_MultiVector> pres_f = Teuchos::null;
      if (pres->HasComponent("face")) {
        pres_f = pres->ViewComponent("face", false);
      }

      int ncells_per = -1;
      for (int col=0; col!=ncols; ++col) {
        auto& col_cells = mesh_->cells_of_column(col);
        if (ncells_per < 0) ncells_per = col_cells.size();
        AMANZI_ASSERT(col_cells.size() == ncells_per);
        auto& col_faces = mesh_->faces_of_column(col);
        AMANZI_ASSERT(col_faces.size() == col_cells.size()+1);
        double z_wt = mesh_->face_centroid(col_faces[0])[z_index] + head_wt;

        if (pres_f.get()) {
          (*pres_f)[0][col_faces[0]] = p_atm + rho * g * head_wt;
          (*flags)[col_faces[0]] = 1;
        }
        for (int lcv_c=0; lcv_c!=col_cells.size(); ++lcv_c) {
          AmanziMesh::Entity_ID c = col_cells[lcv_c];
          AmanziMesh::Entity_ID f = col_faces[lcv_c+1];
          pres_c[0][c] = p_atm + rho * g * (z_wt - mesh_->cell_centroid(c)[z_index]);
          if (pres_f.get()) {
            (*pres_f)[0][f] = p_atm + rho * g * (z_wt - mesh_->face_centroid(f)[z_index]);
            (*flags)[f] = 1;
          }
        }
      }
    }
    if (pres->HasComponent("face")) {
      // communicate, then deal with horizontal-normal faces
      pres->ScatterMasterToGhosted("cell");
      {
        const Epetra_MultiVector& pres_c = *pres->ViewComponent("cell", false);
        Epetra_MultiVector& pres_f = *pres->ViewComponent("face", false);
        for (AmanziMesh::Entity_ID f=0; f!=pres_f.MyLength(); ++f) {
          if (!(*flags)[f]) {
            AmanziMesh::Entity_ID_List f_cells;
            mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &f_cells);
            if (f_cells.size() == 1) {
              // boundary face, use the cell value as the water table is
              // assumed to parallel the cell structure
              pres_f[0][f] = pres_c[0][f_cells[0]];
              (*flags)[f] = 1;
            } else {
              AMANZI_ASSERT(f_cells.size() == 2);
              // interpolate between cells
              pres_f[0][f] = (pres_c[0][f_cells[0]] + pres_c[0][f_cells[1]]) / 2.;
              (*flags)[f] = 1;
            }
          }
        }
      }
    }
    S->GetField(key_, name_)->set_initialized();
  }

  // constant head datum
  if (plist_->sublist("initial condition").isParameter("hydrostatic water level [m]")) {
    double z_wt = plist_->sublist("initial condition").get<double>("hydrostatic water level [m]");
    double rho = plist_->sublist("initial condition").get<double>("hydrostatic water density [kg m^-3]");
    Teuchos::RCP<const Epetra_Vector> gvec = S->GetConstantVectorData("gravity");
    int z_index = mesh_->space_dimension() - 1;
    double g = -(*gvec)[z_index];
    double p_atm = *S->GetScalarData("atmospheric_pressure");

    Teuchos::RCP<CompositeVector> pres = S->GetFieldData(key_, name_);
    Epetra_MultiVector& pres_c = *pres->ViewComponent("cell", false);
    for (int c=0; c!=pres_c.MyLength(); ++c) {
      pres_c[0][c] = p_atm + rho * g * (z_wt - mesh_->cell_centroid(c)[z_index]);
    }

    if (pres->HasComponent("face")) {
      Epetra_MultiVector& pres_f = *pres->ViewComponent("face", false);
      for (int f=0; f!=pres_f.MyLength(); ++f) {
        pres_f[0][f] = p_atm + rho * g * (z_wt - mesh_->face_centroid(f)[z_index]);
      }
    }
    S->GetField(key_, name_)->set_initialized();
  }
}


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void Richards::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  double dt = t_new - t_old;
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Commiting state." << std::endl;

  PK_PhysicalBDF_Default::CommitStep(t_old, t_new, S);

  // update BCs, rel perm
  UpdateBoundaryConditions_(S.ptr());
  bool update = UpdatePermeabilityData_(S.ptr());
  update |= S->GetFieldEvaluator(key_)->HasFieldChanged(S.ptr(), name_);
  update |= S->GetFieldEvaluator(mass_dens_key_)->HasFieldChanged(S.ptr(), name_);

  if (update) {
    // update the stiffness matrix and derive fluxes
    Teuchos::RCP<const CompositeVector> rel_perm = S->GetFieldData(uw_coef_key_);
    Teuchos::RCP<const CompositeVector> rho = S->GetFieldData(mass_dens_key_);
    Teuchos::RCP<CompositeVector> pres = S->GetFieldData(key_, name_);

    matrix_->Init();
    matrix_diff_->SetDensity(rho);
    matrix_diff_->SetScalarCoefficient(rel_perm, Teuchos::null);
    matrix_diff_->UpdateMatrices(Teuchos::null, pres.ptr());
    matrix_diff_->ApplyBCs(true, true, true);

    // derive fluxes
    Teuchos::RCP<CompositeVector> flux = S->GetFieldData(flux_key_, name_);
    matrix_diff_->UpdateFlux(pres.ptr(), flux.ptr());

    if (compute_boundary_values_) {
      applyDirichletBCs(*bc_, *pres);
    }
  }

  // As a diagnostic, calculate the mass balance error
// #if DEBUG_FLAG
//   if (S_next_ != Teuchos::null) {
//     Teuchos::RCP<const CompositeVector> wc1 = S_next_->GetFieldData(conserved_key_);
//     Teuchos::RCP<const CompositeVector> wc0 = S_->GetFieldData(conserved_key_);
//     Teuchos::RCP<const CompositeVector> mass_flux = S->GetFieldData(flux_key_, name_);
//     CompositeVector error(*wc1);

//     for (unsigned int c=0; c!=error.size("cell"); ++c) {
//       error("cell",c) = (*wc1)("cell",c) - (*wc0)("cell",c);

//       AmanziMesh::Entity_ID_List faces;
//       std::vector<int> dirs;
//       mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
//       for (unsigned int n=0; n!=faces.size(); ++n) {
//         error("cell",c) += (*mass_flux)("face",faces[n]) * dirs[n] * dt;
//       }
//     }

//     double einf(0.0);
//     error.NormInf(&einf);

//     // VerboseObject stuff.
//     Teuchos::OSTab tab = vo_->getOSTab();
//     *vo_->os() << "Final Mass Balance Error: " << einf << std::endl;
//   }
// #endif
};


// -----------------------------------------------------------------------------
// Check for controls on saturation
// -----------------------------------------------------------------------------
bool
Richards::ValidStep()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Validating time step." << std::endl;

  if (sat_change_limit_ > 0.0) {
    const Epetra_MultiVector& sl_new = *S_next_->GetFieldData(sat_key_)
        ->ViewComponent("cell",false);
    const Epetra_MultiVector& sl_old = *S_inter_->GetFieldData(sat_key_)
        ->ViewComponent("cell",false);
    Epetra_MultiVector dsl(sl_new);
    dsl.Update(-1., sl_old, 1.);
    double change = 0.;
    dsl.NormInf(&change);

    if (change > sat_change_limit_) {
      if (vo_->os_OK(Teuchos::VERB_LOW))
        *vo_->os() << "Invalid time step, max sl change="
                   << change << " > limit=" << sat_change_limit_ << std::endl;
      return false;
    }
  }
  if (sat_ice_change_limit_ > 0.0) {
    const Epetra_MultiVector& si_new = *S_next_->GetFieldData(sat_ice_key_)
        ->ViewComponent("cell",false);
    const Epetra_MultiVector& si_old = *S_inter_->GetFieldData(sat_ice_key_)
        ->ViewComponent("cell",false);
    Epetra_MultiVector dsi(si_new);
    dsi.Update(-1., si_old, 1.);
    double change = 0.;
    dsi.NormInf(&change);

    if (change > sat_ice_change_limit_) {
      if (vo_->os_OK(Teuchos::VERB_LOW))
        *vo_->os() << "Invalid time step, max si change="
                   << change << " > limit=" << sat_ice_change_limit_ << std::endl;
      return false;
    }
  }
  return PK_PhysicalBDF_Default::ValidStep();
}


// -----------------------------------------------------------------------------
// Update any diagnostic variables prior to vis (in this case velocity field).
// -----------------------------------------------------------------------------
void Richards::CalculateDiagnostics(const Teuchos::RCP<State>& S)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Calculating diagnostic variables." << std::endl;

  // update the cell velocities
  UpdateBoundaryConditions_(S.ptr());

  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(key_);
  Teuchos::RCP<const CompositeVector> rel_perm = S->GetFieldData(uw_coef_key_);
  Teuchos::RCP<const CompositeVector> rho = S->GetFieldData(mass_dens_key_);
  // update the stiffness matrix
  matrix_diff_->SetDensity(rho);
  matrix_diff_->SetScalarCoefficient(rel_perm, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, pres.ptr());
  matrix_diff_->ApplyBCs(true, true, true);

  // derive fluxes
  Teuchos::RCP<CompositeVector> flux = S->GetFieldData(flux_key_, name_);
  matrix_diff_->UpdateFlux(pres.ptr(), flux.ptr());

  UpdateVelocity_(S.ptr());
};


// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
bool Richards::UpdatePermeabilityData_(const Teuchos::Ptr<State>& S)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating permeability?";

  Teuchos::RCP<const CompositeVector> rel_perm = S->GetFieldData(coef_key_);

  const Epetra_MultiVector& rel_perm_c =
      *S->GetFieldData(coef_key_)->ViewComponent("cell",false);

  const Epetra_MultiVector& pre_c =
      *S->GetFieldData(key_)->ViewComponent("cell",false);    

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

   for (int c = 0; c < ncells_owned; c++) {
     std::cout << "rel_perm_c[0][" << c << "] = " << rel_perm_c[0][c] << std::endl;
     std::cout << "pre_c[0][" << c << "] = " << pre_c[0][c] << std::endl;
     pre_c[0][c] = 101325.;
   }

  std::cout << "coef_key_ = " << coef_key_ << std::endl;
  bool update_perm = S->GetFieldEvaluator(coef_key_)
      ->HasFieldChanged(S, name_);

  // requirements due to the upwinding method
  if (Krel_method_ == Operators::UPWIND_METHOD_TOTAL_FLUX) {
    bool update_dir = S->GetFieldEvaluator(mass_dens_key_)
        ->HasFieldChanged(S, name_);
    update_dir |= S->GetFieldEvaluator(key_)->HasFieldChanged(S, name_);

    if (update_dir) {
      // update the direction of the flux -- note this is NOT the flux
      Teuchos::RCP<const CompositeVector> rho = S->GetFieldData(mass_dens_key_);
      Teuchos::RCP<CompositeVector> flux_dir = S->GetFieldData(flux_dir_key_, name_);
      Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(key_);


      face_matrix_diff_->SetDensity(rho);
      face_matrix_diff_->UpdateMatrices(Teuchos::null, pres.ptr());
      //if (!pres->HasComponent("face"))
      face_matrix_diff_->ApplyBCs(true, true, true);
      face_matrix_diff_->UpdateFlux(pres.ptr(), flux_dir.ptr());

      if (clobber_boundary_flux_dir_) {
        Epetra_MultiVector& flux_dir_f = *flux_dir->ViewComponent("face",false);

        auto& markers = bc_markers();
        auto& values = bc_values();

        for (int f=0; f!=markers.size(); ++f) {
          if (markers[f] == Operators::OPERATOR_BC_NEUMANN) {
            AmanziMesh::Entity_ID_List cells;
            mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
            AMANZI_ASSERT(cells.size() == 1);
            int c = cells[0];
            AmanziMesh::Entity_ID_List faces;
            std::vector<int> dirs;
            mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
            int i = std::find(faces.begin(), faces.end(), f) - faces.begin();

            flux_dir_f[0][f] = values[f]*dirs[i];
          }
        }
      }
    }

    update_perm |= update_dir;
  }

  if (update_perm) {
    Teuchos::RCP<CompositeVector> uw_rel_perm = S->GetFieldData(uw_coef_key_, name_);

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

    if (clobber_policy_ == "clobber") {
      Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->ViewComponent("face",false);
      uw_rel_perm_f.Export(rel_perm_bf, vandelay, Insert);
    } else if (clobber_policy_ == "max") {
      Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->ViewComponent("face",false);
      const auto& fmap = mesh_->face_map(true);
      const auto& bfmap = mesh_->exterior_face_map(true);
      for (int bf=0; bf!=rel_perm_bf.MyLength(); ++bf) {
        auto f = fmap.LID(bfmap.GID(bf));
        if (rel_perm_bf[0][bf] > uw_rel_perm_f[0][f]) {
          uw_rel_perm_f[0][f] = rel_perm_bf[0][bf];
        }
      }
    } else if (clobber_policy_ == "unsaturated") {
      // clobber only when the interior cell is unsaturated
      Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->ViewComponent("face",false);
      const Epetra_MultiVector& pres = *S->GetFieldData(key_)->ViewComponent("cell",false);
      const auto& fmap = mesh_->face_map(true);
      const auto& bfmap = mesh_->exterior_face_map(true);
      for (int bf=0; bf!=rel_perm_bf.MyLength(); ++bf) {
        auto f = fmap.LID(bfmap.GID(bf));
        AmanziMesh::Entity_ID_List fcells;
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &fcells);
        AMANZI_ASSERT(fcells.size() == 1);
        if (pres[0][fcells[0]] < 101225.) {
          uw_rel_perm_f[0][f] = rel_perm_bf[0][bf];
        } else if (pres[0][fcells[0]] < 101325.) {
          double frac = (101325. - pres[0][fcells[0]])/100.;
          uw_rel_perm_f[0][f] = rel_perm_bf[0][bf] * frac + uw_rel_perm_f[0][f] * (1-frac);
        }
      }
    }

    if (uw_rel_perm->HasComponent("face"))
      uw_rel_perm->ScatterMasterToGhosted("face");
  }

  // debugging
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << " " << update_perm << std::endl;
  }
  return update_perm;
};


bool Richards::UpdatePermeabilityDerivativeData_(const Teuchos::Ptr<State>& S)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating permeability derivatives?";

  bool update_perm = S->GetFieldEvaluator(coef_key_)->HasFieldDerivativeChanged(S, name_, key_);
  Teuchos::RCP<const CompositeVector> drel_perm = S->GetFieldData(dcoef_key_);

  if (update_perm) {
    if (!duw_coef_key_.empty()) {
      Teuchos::RCP<CompositeVector> duw_rel_perm = S->GetFieldData(duw_coef_key_, name_);
      duw_rel_perm->PutScalar(0.);

      // Upwind, only overwriting boundary faces if the wind says to do so.
      upwinding_deriv_->Update(S);

      duw_rel_perm->ScatterMasterToGhosted("face");
    } else {
      drel_perm->ScatterMasterToGhosted("cell");
    }
  }

  // debugging
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
    *vo_->os() << " " << update_perm << std::endl;
  }
  return update_perm;
};


// -----------------------------------------------------------------------------
// Compute boundary condition functions at the current time.
// -----------------------------------------------------------------------------
void Richards::ComputeBoundaryConditions_(const Teuchos::Ptr<State>& S)
{
  bc_pressure_->Compute(S->time());
  bc_head_->Compute(S->time());
  bc_level_->Compute(S->time());
  bc_flux_->Compute(S->time());
  bc_seepage_->Compute(S->time());
  bc_seepage_infilt_->Compute(S->time());
}


// -----------------------------------------------------------------------------
// Push boundary conditions into the global array.
// -----------------------------------------------------------------------------
void Richards::UpdateBoundaryConditions_(const Teuchos::Ptr<State>& S, bool kr)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Updating BCs." << std::endl;

  auto& markers = bc_markers();
  auto& values = bc_values();

  // initialize all to 0
  for (unsigned int n=0; n!=markers.size(); ++n) {
    markers[n] = Operators::OPERATOR_BC_NONE;
    values[n] = 0.0;
  }

  // count for debugging
  std::vector<int> bc_counts;
  std::vector<std::string> bc_names;

  // Dirichlet-type boundary conditions
  // -------------------------------------
  // pressure boundary conditions -- the primary
  bc_counts.push_back(bc_pressure_->size());
  bc_names.push_back(key_);
  for (const auto& bc : *bc_pressure_) {
    int f = bc.first;
#ifdef ENABLE_DBC
    AmanziMesh::Entity_ID_List cells;
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    AMANZI_ASSERT(cells.size() == 1);
#endif
    markers[f] = Operators::OPERATOR_BC_DIRICHLET;
    values[f] = bc.second;
  }

  // head boundary conditions
  bc_counts.push_back(bc_head_->size());
  bc_names.push_back("head");
  if (bc_head_->size() > 0) {
    double p_atm = *S->GetScalarData("atmospheric_pressure");
    int z_index = mesh_->space_dimension() - 1;
    double g = -(*S->GetConstantVectorData("gravity"))[z_index];

    for (const auto& bc : *bc_head_) {
      int f = bc.first;

      AmanziMesh::Entity_ID_List cells;
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      AMANZI_ASSERT(cells.size() == 1);

      // we need to find the elevation of the surface, but finding the top edge
      // of this stack of faces is not possible currently.  The best approach
      // is instead to work with the cell.
      int col = mesh_->column_ID(cells[0]);
      double z_surf = mesh_->face_centroid(mesh_->faces_of_column(col)[0])[z_index];
      double z_wt = bc.second + z_surf;

      markers[f] = Operators::OPERATOR_BC_DIRICHLET;
      // note, here the cell centroid's z is used to relate to the column's top
      // face centroid, specifically NOT the boundary face's centroid.
      values[f] = p_atm + bc_rho_water_ * g *
        (z_wt - mesh_->cell_centroid(cells[0])[z_index]);
    }
  }

  // fixed level head boundary conditions
  bc_counts.push_back(bc_head_->size());
  bc_names.push_back("head");
  if (bc_level_->size() > 0) {
    double p_atm = *S->GetScalarData("atmospheric_pressure");
    int z_index = mesh_->space_dimension() - 1;
    double g = -(*S->GetConstantVectorData("gravity"))[z_index];

    for (const auto& bc : *bc_level_) {
      int f = bc.first;
      markers[f] = Operators::OPERATOR_BC_DIRICHLET;
      values[f] = p_atm + bc_rho_water_ * g *
        (bc.second - mesh_->face_centroid(f)[z_index]);
    }
  }

  // Neumann type boundary conditions
  // -------------------------------------
  const Epetra_MultiVector& rel_perm =
    *S->GetFieldData(uw_coef_key_)->ViewComponent("face",false);

  // standard Neumann flux BCs
  bc_counts.push_back(bc_flux_->size());
  bc_names.push_back(flux_key_);

  if (!infiltrate_only_if_unfrozen_) {
    for (const auto& bc : *bc_flux_) {
      int f = bc.first;
#ifdef ENABLE_DBC
    AmanziMesh::Entity_ID_List cells;
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    AMANZI_ASSERT(cells.size() == 1);
#endif
      markers[f] = Operators::OPERATOR_BC_NEUMANN;
      values[f] = bc.second;
      if (!kr && rel_perm[0][f] > 0.) values[f] /= rel_perm[0][f];
    }

  } else {
    // Neumann boundary conditions that turn off if temp < freezing
    const Epetra_MultiVector& temp = *S->GetFieldData(Keys::getKey(domain_,"temperature"))->ViewComponent("face");
    for (const auto& bc : *bc_flux_) {
      int f = bc.first;
#ifdef ENABLE_DBC
    AmanziMesh::Entity_ID_List cells;
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    AMANZI_ASSERT(cells.size() == 1);
#endif
      markers[f] = Operators::OPERATOR_BC_NEUMANN;
      if (temp[0][f] > 273.15) {
        values[f] = bc.second;
        if (!kr && rel_perm[0][f] > 0.) values[f] /= rel_perm[0][f];
      } else {
        values[f] = 0.;
      }
    }
  }

  // seepage face -- pressure <= specified value (usually 101325), outward mass flux >= 0
  S->GetFieldData(flux_key_)->ScatterMasterToGhosted("face");
  const Epetra_MultiVector& flux = *S->GetFieldData(flux_key_)->ViewComponent("face", true);

  const double& p_atm = *S->GetScalarData("atmospheric_pressure");
  Teuchos::RCP<const CompositeVector> u = S->GetFieldData(key_);
  double seepage_tol = 10.;

  bc_counts.push_back(bc_seepage_->size());
  bc_names.push_back("standard seepage");
  for (const auto& bc : *bc_seepage_) {
    int f = bc.first;
#ifdef ENABLE_DBC
    AmanziMesh::Entity_ID_List cells;
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    AMANZI_ASSERT(cells.size() == 1);
#endif

    //double boundary_pressure = std::max(getFaceOnBoundaryValue(f, *u, *bc_), 101325.); // does not make sense to seep from nonsaturated cells
    double boundary_pressure = getFaceOnBoundaryValue(f, *u, *bc_); // does not make sense to seep from nonsaturated cells
    double boundary_flux = flux[0][f]*getBoundaryDirection(*mesh_, f);
    if (boundary_pressure > bc.second) {
      markers[f] = Operators::OPERATOR_BC_DIRICHLET;
      values[f] = bc.second;
    } else if (boundary_pressure < bc.second - seepage_tol) {
      markers[f] = Operators::OPERATOR_BC_NEUMANN;
      values[f] = 0.;
    } else if (boundary_flux >= 0.) {
      markers[f] = Operators::OPERATOR_BC_DIRICHLET;
      values[f] = bc.second;
    } else {
      markers[f] = Operators::OPERATOR_BC_NEUMANN;
      values[f] = 0.;
    }
  }

  // seepage face -- pressure <= p_atm, outward mass flux is specified
  bc_counts.push_back(bc_seepage_infilt_->size());
  bc_names.push_back("seepage with infiltration");
  {
    Teuchos::Ptr<State> Sl;
    if (bc_seepage_infilt_explicit_ && S_inter_.get()) {
      Sl = S_inter_.ptr();
    } else {
      Sl = S;
    }
    const Epetra_MultiVector& flux = *Sl->GetFieldData(flux_key_)->ViewComponent("face", true);
    Teuchos::RCP<const CompositeVector> u = Sl->GetFieldData(key_);
  int i = 0;
  for (const auto& bc : *bc_seepage_infilt_) {
    int f = bc.first;
#ifdef ENABLE_DBC
    AmanziMesh::Entity_ID_List cells;
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    AMANZI_ASSERT(cells.size() == 1);
#endif

    double flux_seepage_tol = std::abs(bc.second) * .001;
    double boundary_pressure = getFaceOnBoundaryValue(f, *u, *bc_);
    double boundary_flux = flux[0][f]*getBoundaryDirection(*mesh_, f);

    if (i == 0)
      std::cout << "BFlux = " << boundary_flux << " with constraint = " << bc.second - flux_seepage_tol << std::endl;

    if (boundary_flux < bc.second - flux_seepage_tol &&
        boundary_pressure > p_atm + seepage_tol) {
      // both constraints are violated, either option should push things in the right direction
      markers[f] = Operators::OPERATOR_BC_DIRICHLET;
      values[f] = p_atm;
      if (i == 0)
        std::cout << "BC PRESSURE ON SEEPAGE = " << boundary_pressure << " with flux " << flux[0][f]*getBoundaryDirection(*mesh_, f) << " resulted in DIRICHLET pressure " << p_atm << std::endl;

    } else if (boundary_flux >= bc.second - flux_seepage_tol &&
        boundary_pressure > p_atm - seepage_tol) {
      // max pressure condition violated
      markers[f] = Operators::OPERATOR_BC_DIRICHLET;
      values[f] = p_atm;
    if (i == 0)
      std::cout << "BC PRESSURE ON SEEPAGE = " << boundary_pressure << " with flux " << flux[0][f]*getBoundaryDirection(*mesh_, f) << " resulted in DIRICHLET pressure " << p_atm << std::endl;

    } else if (boundary_flux < bc.second - flux_seepage_tol &&
        boundary_pressure <= p_atm + seepage_tol) {
      // max infiltration violated
      markers[f] = Operators::OPERATOR_BC_NEUMANN;
      values[f] = bc.second;
    if (i == 0)
      std::cout << "BC PRESSURE ON SEEPAGE = " << boundary_pressure << " with flux " << flux[0][f]*getBoundaryDirection(*mesh_, f) << " resulted in NEUMANN flux " << bc.second << std::endl;

    } else if (boundary_flux >= bc.second - flux_seepage_tol &&
        boundary_pressure <= p_atm - seepage_tol) {
      // both conditions are valid
      markers[f] = Operators::OPERATOR_BC_NEUMANN;
      values[f] = bc.second;
    if (i == 0)
      std::cout << "BC PRESSURE ON SEEPAGE = " << boundary_pressure << " with flux " << flux[0][f]*getBoundaryDirection(*mesh_, f) << " resulted in NEUMANN flux " << bc.second << std::endl;

    } else {
      AMANZI_ASSERT(0);
    }
    i++;
  }
  }

  // surface coupling
  bc_counts.push_back(0);
  bc_names.push_back("surface coupling (head)");

  if (coupled_to_surface_via_head_) {
    // Face is Dirichlet with value of surface head
    Teuchos::RCP<const AmanziMesh::Mesh> surface = S->GetMesh("surface");
    const Epetra_MultiVector& head = *S->GetFieldData("surface_pressure")
        ->ViewComponent("cell",false);

    unsigned int ncells_surface = head.MyLength();
    bc_counts[bc_counts.size()-1] = ncells_surface;

    for (unsigned int c=0; c!=ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face
      AmanziMesh::Entity_ID f = surface->entity_get_parent(AmanziMesh::CELL, c);
#ifdef ENABLE_DBC
      AmanziMesh::Entity_ID_List cells;
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      AMANZI_ASSERT(cells.size() == 1);
#endif
      // -- set that value to dirichlet
      markers[f] = Operators::OPERATOR_BC_DIRICHLET;
      values[f] = head[0][c];
    }
  }

  // surface coupling
  bc_counts.push_back(0);
  bc_names.push_back("surface coupling (flux)");
  if (coupled_to_surface_via_flux_) {
    // Face is Neumann with value of surface residual
    Teuchos::RCP<const AmanziMesh::Mesh> surface = S->GetMesh(Keys::getDomain(ss_flux_key_));
    const Epetra_MultiVector& ss_flux = *S->GetFieldData(ss_flux_key_)->ViewComponent("cell",false);
    unsigned int ncells_surface = ss_flux.MyLength();
    bc_counts[bc_counts.size()-1] = ncells_surface;
    for (unsigned int c=0; c!=ncells_surface; ++c) {
      // -- get the surface cell's equivalent subsurface face
      AmanziMesh::Entity_ID f = surface->entity_get_parent(AmanziMesh::CELL, c);
#ifdef ENABLE_DBC
      AmanziMesh::Entity_ID_List cells;
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      AMANZI_ASSERT(cells.size() == 1);
#endif
      // -- set that value to Neumann
      markers[f] = Operators::OPERATOR_BC_NEUMANN;

      // NOTE: the flux provided by the coupler is in units of mols / s, where
      //       as Neumann BCs are in units of mols / s / A.  The right A must
      //       be chosen, as it is the subsurface mesh's face area, not the
      //       surface mesh's cell area.
      values[f] = ss_flux[0][c] / mesh_->face_area(f);

      if (!kr && rel_perm[0][f] > 0.) values[f] /= rel_perm[0][f];
    }
  }

  // mark all remaining boundary conditions as zero flux conditions
  AmanziMesh::Entity_ID_List cells;
  int n_default = 0;
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    if (markers[f] == Operators::OPERATOR_BC_NONE) {
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int ncells = cells.size();

      if (ncells == 1) {
        n_default++;
        markers[f] = Operators::OPERATOR_BC_NEUMANN;
        values[f] = 0.0;
      }
    }
  }
  bc_names.push_back("default (zero flux)");
  bc_counts.push_back(n_default);

  // report on counts
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    std::vector<int> bc_counts_global(bc_counts.size(), 0);
    mesh_->get_comm()->SumAll(&bc_counts[0], &bc_counts_global[0], bc_counts.size());

    *vo_->os() << "  BCs applied:" << std::endl;

    for (int i=0; i!=bc_counts_global.size(); ++i) {
      *vo_->os() << "    " << bc_names[i] << ": " << bc_counts_global[i] << std::endl;
    }
  }

};


bool Richards::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
        Teuchos::RCP<TreeVector> u)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Modifying predictor:" << std::endl;

  // update boundary conditions
  ComputeBoundaryConditions_(S_next_.ptr());
  UpdateBoundaryConditions_(S_next_.ptr());
  db_->WriteBoundaryConditions(bc_markers(), bc_values());

  // push Dirichlet data into predictor
  applyDirichletBCs(*bc_, *u->Data());

  bool changed(false);
  if (modify_predictor_bc_flux_ ||
      (modify_predictor_first_bc_flux_ &&
       ((S_next_->cycle() == 0) || (S_next_->cycle() == 1)))) {
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

bool Richards::ModifyPredictorFluxBCs_(double h, Teuchos::RCP<TreeVector> u)
{
  if (!u->Data()->HasComponent("face")) return false;

  auto& markers = bc_markers();
  auto& values = bc_values();

  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  modifications to deal with nonlinearity at flux BCs" << std::endl;

  if (flux_predictor_ == Teuchos::null) {
    flux_predictor_ = Teuchos::rcp(new PredictorDelegateBCFlux(S_next_, mesh_, matrix_diff_,
            wrms_, &markers, &values));
  }

  UpdatePermeabilityData_(S_next_.ptr());
  Teuchos::RCP<const CompositeVector> rel_perm =
    S_next_->GetFieldData(uw_coef_key_);

  matrix_->Init();
  matrix_diff_->SetScalarCoefficient(rel_perm, Teuchos::null);
  Teuchos::RCP<const CompositeVector> rho = S_next_->GetFieldData(mass_dens_key_);
  Teuchos::RCP<const CompositeVector> pres = S_next_->GetFieldData(key_);
  matrix_diff_->SetDensity(rho);
  matrix_diff_->UpdateMatrices(Teuchos::null, pres.ptr());
  //matrix_diff_->ApplyBCs(true, true, true);

  flux_predictor_->ModifyPredictor(h, u);
  ChangedSolution(); // mark the solution as changed, as modifying with
                      // consistent faces will then get the updated boundary
                      // conditions
  return true;
}

bool Richards::ModifyPredictorConsistentFaces_(double h, Teuchos::RCP<TreeVector> u)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  modifications for consistent face pressures." << std::endl;

  CalculateConsistentFaces(u->Data().ptr());
  return true;
}

bool Richards::ModifyPredictorWC_(double h, Teuchos::RCP<TreeVector> u)
{
  AMANZI_ASSERT(0);
  return false;
}


// void Richards::CalculateConsistentFacesForInfiltration_(
//     const Teuchos::Ptr<CompositeVector>& u) {

//  auto& markers = bc_markers();
//  auto& values = bc_values();

//   if (vo_->os_OK(Teuchos::VERB_EXTREME))
//     *vo_->os() << "  modifications to deal with nonlinearity at flux BCs" << std::endl;

//   if (flux_predictor_ == Teuchos::null) {
//     flux_predictor_ = Teuchos::rcp(new PredictorDelegateBCFlux(S_next_, mesh_, matrix_,
//             wrms_, &markers, &values));
//   }

//   // update boundary conditions
//   bc_pressure_->Compute(S_next_->time());
//   bc_flux_->Compute(S_next_->time());
//   UpdateBoundaryConditions_(S_next_.ptr());

//   bool update = UpdatePermeabilityData_(S_next_.ptr());
//   Teuchos::RCP<const CompositeVector> rel_perm =
//       S_next_->GetFieldData(uw_coef_key_);
//   matrix_->CreateMFDstiffnessMatrices(rel_perm.ptr());
//   matrix_->CreateMFDrhsVectors();
//   Teuchos::RCP<const CompositeVector> rho = S_next_->GetFieldData(mass_dens_key_);
//   Teuchos::RCP<const Epetra_Vector> gvec = S_next_->GetConstantVectorData("gravity");
//   AddGravityFluxes_(gvec.ptr(), rel_perm.ptr(), rho.ptr(), matrix_.ptr());
//   matrix_->ApplyBoundaryConditions(markers, values);

//   flux_predictor_->ModifyPredictor(u);
// }

void Richards::CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u)
{
  if (!u->HasComponent("face")) return; // not need

  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Modifying predictor for consistent faces" << std::endl;


  // average cells to faces to give a reasonable initial guess

  u->ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& u_c = *u->ViewComponent("cell",true);
  Epetra_MultiVector& u_f = *u->ViewComponent("face",false);

  int f_owned = u_f.MyLength();
  for (int f=0; f!=f_owned; ++f) {
    AmanziMesh::Entity_ID_List cells;
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    double face_value = 0.0;
    for (int n=0; n!=ncells; ++n) {
      face_value += u_c[0][cells[n]];
    }
    u_f[0][f] = face_value / ncells;
  }
  ChangedSolution();

  // Using the old BCs, so should use the old rel perm?
  // update the rel perm according to the scheme of choice
  //  UpdatePermeabilityData_(S_next_.ptr());

  Teuchos::RCP<const CompositeVector> rel_perm =
    S_next_->GetFieldData(uw_coef_key_);

  S_next_->GetFieldEvaluator(mass_dens_key_)
      ->HasFieldChanged(S_next_.ptr(), name_);
  Teuchos::RCP<const CompositeVector> rho =
      S_next_->GetFieldData(mass_dens_key_);


  // Update the preconditioner with darcy and gravity fluxes
  matrix_->Init();
  matrix_diff_->SetDensity(rho);
  matrix_diff_->SetScalarCoefficient(rel_perm, Teuchos::null);
  matrix_diff_->UpdateMatrices(Teuchos::null, u);
  matrix_diff_->ApplyBCs(true, true, true);

  // derive the consistent faces, involves a solve

  db_->WriteVector(" p_cf guess:", u.ptr(), true);
  matrix_diff_->UpdateConsistentFaces(*u);
  db_->WriteVector(" p_cf soln:", u.ptr(), true);

}


/* ******************************************************************
* Clip pressure using pressure threshold.
****************************************************************** */
void Richards::ClipHydrostaticPressure(double pmin, Epetra_MultiVector& p)
{
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int c = 0; c < ncells_owned; c++) p[0][c] = std::max(p[0][c], pmin);
}


// -----------------------------------------------------------------------------
// Check admissibility of the solution guess.
// -----------------------------------------------------------------------------
bool Richards::IsAdmissible(Teuchos::RCP<const TreeVector> up)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "  Checking admissibility..." << std::endl;

  // For some reason, wandering PKs break most frequently with an unreasonable
  // pressure.  This simply tries to catch that before it happens.
  Teuchos::RCP<const CompositeVector> pres = up->Data();
  double minT, maxT;

  const Epetra_MultiVector& pres_c = *pres->ViewComponent("cell",false);
  double minT_c(1.e15), maxT_c(-1.e15);
  int min_c(-1), max_c(-1);
  for (int c=0; c!=pres_c.MyLength(); ++c) {
    if (pres_c[0][c] < minT_c) {
      minT_c = pres_c[0][c];
      min_c = c;
    }
    if (pres_c[0][c] > maxT_c) {
      maxT_c = pres_c[0][c];
      max_c = c;
    }
  }

  double minT_f(1.e15), maxT_f(-1.e15);
  int min_f(-1), max_f(-1);
  if (pres->HasComponent("face")) {
    const Epetra_MultiVector& pres_f = *pres->ViewComponent("face",false);
    for (int f=0; f!=pres_f.MyLength(); ++f) {
      if (pres_f[0][f] < minT_f) {
        minT_f = pres_f[0][f];
        min_f = f;
      }
      if (pres_f[0][f] > maxT_f) {
        maxT_f = pres_f[0][f];
        max_f = f;
      }
    }
    minT = std::min(minT_c, minT_f);
    maxT = std::max(maxT_c, maxT_f);

  } else {
    minT = minT_c;
    maxT = maxT_c;
  }

  double minT_l = minT;
  double maxT_l = maxT;
  mesh_->get_comm()->MaxAll(&maxT_l, &maxT, 1);
  mesh_->get_comm()->MinAll(&minT_l, &minT, 1);

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "    Admissible p? (min/max): " << minT << ",  " << maxT << std::endl;
  }


  Teuchos::RCP<const Comm_type> comm_p = mesh_->get_comm();
    Teuchos::RCP<const MpiComm_type> mpi_comm_p =
      Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm_p);
    const MPI_Comm& comm = mpi_comm_p->Comm();

  if (minT < -1.e9 || maxT > 1.e8) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << " is not admissible, as it is not within bounds of constitutive models:" << std::endl;
      ENorm_t global_minT_c, local_minT_c;
      ENorm_t global_maxT_c, local_maxT_c;

      local_minT_c.value = minT_c;
      local_minT_c.gid = pres_c.Map().GID(min_c);
      local_maxT_c.value = maxT_c;
      local_maxT_c.gid = pres_c.Map().GID(max_c);

      MPI_Allreduce(&local_minT_c, &global_minT_c, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
      MPI_Allreduce(&local_maxT_c, &global_maxT_c, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
      *vo_->os() << "   cells (min/max): [" << global_minT_c.gid << "] " << global_minT_c.value
                 << ", [" << global_maxT_c.gid << "] " << global_maxT_c.value << std::endl;

      if (pres->HasComponent("face")) {
        const Epetra_MultiVector& pres_f = *pres->ViewComponent("face",false);
        ENorm_t global_minT_f, local_minT_f;
        ENorm_t global_maxT_f, local_maxT_f;

        local_minT_f.value = minT_f;
        local_minT_f.gid = pres_f.Map().GID(min_f);
        local_maxT_f.value = maxT_f;
        local_maxT_f.gid = pres_f.Map().GID(max_f);


        MPI_Allreduce(&local_minT_f, &global_minT_f, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
        MPI_Allreduce(&local_maxT_f, &global_maxT_f, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
        *vo_->os() << "   cells (min/max): [" << global_minT_f.gid << "] " << global_minT_f.value;
        MPI_Allreduce(&local_minT_f, &global_minT_f, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
        MPI_Allreduce(&local_maxT_f, &global_maxT_f, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        *vo_->os() << "   faces (min/max): [" << global_minT_f.gid << "] " << global_minT_f.value
                   << ", [" << global_maxT_f.gid << "] " << global_maxT_f.value << std::endl;
      }
    }
    return false;
  }
  return true;
}


AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
Richards::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                 Teuchos::RCP<const TreeVector> u,
                 Teuchos::RCP<TreeVector> du)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // if the primary variable has boundary face, this is for upwinding rel
  // perms and is never actually used.  Make sure it does not go to undefined
  // pressures.
  if (du->Data()->HasComponent("boundary_face"))  {
    du->Data()->ViewComponent("boundary_face")->PutScalar(0.);
  }

  // debugging -- remove me! --etc
  for (CompositeVector::name_iterator comp=du->Data()->begin();
       comp!=du->Data()->end(); ++comp) {
    Epetra_MultiVector& du_c = *du->Data()->ViewComponent(*comp,false);
    double max, l2;
    du_c.NormInf(&max);
    du_c.Norm2(&l2);
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "Linf, L2 pressure correction (" << *comp << ") = " << max << ", " << l2 << std::endl;
    }
  }

  // limit by capping corrections when they cross atmospheric pressure
  // (where pressure derivatives are discontinuous)
  int my_limited = 0;
  int n_limited_spurt = 0;
  if (patm_limit_ > 0.) {
    double patm = *S_next_->GetScalarData("atmospheric_pressure");
    for (CompositeVector::name_iterator comp=du->Data()->begin();
         comp!=du->Data()->end(); ++comp) {
      Epetra_MultiVector& du_c = *du->Data()->ViewComponent(*comp,false);
      const Epetra_MultiVector& u_c = *u->Data()->ViewComponent(*comp,false);

      for (int c=0; c!=du_c.MyLength(); ++c) {
        if ((u_c[0][c] < patm) &&
            (u_c[0][c] - du_c[0][c] > patm + patm_limit_)) {
          du_c[0][c] = u_c[0][c] - (patm + patm_limit_);
          my_limited++;
        } else if ((u_c[0][c] > patm) &&
                   (u_c[0][c] - du_c[0][c] < patm - patm_limit_)) {
          du_c[0][c] = u_c[0][c] - (patm - patm_limit_);
          my_limited++;
        }
      }
    }
    mesh_->get_comm()->MaxAll(&my_limited, &n_limited_spurt, 1);
  }

  if (n_limited_spurt > 0) {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "  limiting the spurt." << std::endl;
    }
  }

  // debugging -- remove me! --etc
  for (CompositeVector::name_iterator comp=du->Data()->begin();
       comp!=du->Data()->end(); ++comp) {
    Epetra_MultiVector& du_c = *du->Data()->ViewComponent(*comp,false);
    double max, l2;
    du_c.NormInf(&max);
    du_c.Norm2(&l2);
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "Linf, L2 pressure correction (" << *comp << ") = " << max << ", " << l2 << std::endl;
    }
  }

  // Limit based on a max pressure change
  my_limited = 0;
  int n_limited_change = 0;
  if (p_limit_ >= 0.) {
    for (CompositeVector::name_iterator comp=du->Data()->begin();
         comp!=du->Data()->end(); ++comp) {
      Epetra_MultiVector& du_c = *du->Data()->ViewComponent(*comp,false);

      double max;
      du_c.NormInf(&max);
      if (vo_->os_OK(Teuchos::VERB_HIGH)) {
        *vo_->os() << "Max pressure correction (" << *comp << ") = " << max << std::endl;
      }

      for (int c=0; c!=du_c.MyLength(); ++c) {
        if (std::abs(du_c[0][c]) > p_limit_) {
          du_c[0][c] = ((du_c[0][c] > 0) - (du_c[0][c] < 0)) * p_limit_;
          my_limited++;
        }
      }
    }

    mesh_->get_comm()->MaxAll(&my_limited, &n_limited_change, 1);
  }

  // debugging -- remove me! --etc
  for (CompositeVector::name_iterator comp=du->Data()->begin();
       comp!=du->Data()->end(); ++comp) {
    Epetra_MultiVector& du_c = *du->Data()->ViewComponent(*comp,false);
    double max, l2;
    du_c.NormInf(&max);
    du_c.Norm2(&l2);
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "Linf, L2 pressure correction (" << *comp << ") = " << max << ", " << l2 << std::endl;
    }
  }

  if (n_limited_change > 0) {
    if (vo_->os_OK(Teuchos::VERB_HIGH)) {
      *vo_->os() << "  limited by pressure." << std::endl;
    }
  }

  if (n_limited_spurt > 0) {
    return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED_LAG_BACKTRACKING;
  } else if (n_limited_change > 0) {
    return AmanziSolvers::FnBaseDefs::CORRECTION_MODIFIED;
  }

  return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
}

/* ******************************************************************
* Return a pointer to a local operator
****************************************************************** */
Teuchos::RCP<Operators::Operator> Richards::my_operator(
    const Operators::OperatorType& type)
{
  if (type == Operators::OPERATOR_MATRIX) return matrix_;
  else if (type == Operators::OPERATOR_PRECONDITIONER_RAW) return preconditioner_;
  return Teuchos::null;
}

} // namespace
} // namespace
