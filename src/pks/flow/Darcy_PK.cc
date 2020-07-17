/*
  Flow PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

// TPLs
#include "Epetra_Import.h"
#include "Epetra_Vector.h"

// Amanzi
#include "constant_variable_field_evaluator.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "LinearOperatorFactory.hh"
#include "PDE_DiffusionFactory.hh"
#include "PDE_DiffusionFracturedMatrix.hh"
#include "primary_variable_field_evaluator.hh"
#include "TimestepControllerFactory.hh"
#include "Tensor.hh"
#include "UniqueLocalIndex.hh"
#include "WhetStoneMeshUtils.hh"

// Amanzi::Flow
#include "Darcy_PK.hh"
#include "DarcyVelocityEvaluator.hh"
#include "FlowDefs.hh"
#include "FracturePermModelPartition.hh"
#include "FracturePermModelEvaluator.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* New constructor: extracts lists and requires fields.
****************************************************************** */
Darcy_PK::Darcy_PK(Teuchos::ParameterList& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& glist,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& soln) :
  Flow_PK(pk_tree, glist, S, soln),
  soln_(soln)
{
  S_ = S;

  std::string pk_name = pk_tree.name();
  auto found = pk_name.rfind("->");
  if (found != std::string::npos) pk_name.erase(0, found + 2);

  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  fp_list_ = Teuchos::sublist(pk_list, pk_name, true);

  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(fp_list_, "time integrator", true);

  // computational domain
  domain_ = fp_list_->template get<std::string>("domain name", "domain");
}


/* ******************************************************************
* Old constructor for unit tests.
****************************************************************** */
Darcy_PK::Darcy_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                   const std::string& pk_list_name,
                   Teuchos::RCP<State> S,
                   const Teuchos::RCP<TreeVector>& soln) :
    Flow_PK(),
    soln_(soln)
{
  S_ = S;

  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  fp_list_ = Teuchos::sublist(pk_list, pk_list_name, true);

  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(fp_list_, "time integrator", true);

  // domain name
  domain_ = fp_list_->template get<std::string>("domain name", "domain");
}


/* ******************************************************************
* Clean memory.
****************************************************************** */
Darcy_PK::~Darcy_PK()
{
  if (vo_ != Teuchos::null) vo_ = Teuchos::null;
}


/* ******************************************************************
* Define structure of this PK.
****************************************************************** */
void Darcy_PK::Setup(const Teuchos::Ptr<State>& S)
{
  dt_ = -1.0;
  mesh_ = S->GetMesh(domain_);
  dim = mesh_->space_dimension();

  // generate keys here to be available for setup of the base class
  pressure_key_ = Keys::getKey(domain_, "pressure"); 
  hydraulic_head_key_ = Keys::getKey(domain_, "hydraulic_head"); 

  darcy_flux_key_ = Keys::getKey(domain_, "darcy_flux"); 
  darcy_velocity_key_ = Keys::getKey(domain_, "darcy_velocity"); 

  permeability_key_ = Keys::getKey(domain_, "permeability"); 
  porosity_key_ = Keys::getKey(domain_, "porosity"); 

  specific_yield_key_ = Keys::getKey(domain_, "specific_yield"); 
  specific_storage_key_ = Keys::getKey(domain_, "specific_storage"); 
  saturation_liquid_key_ = Keys::getKey(domain_, "saturation_liquid"); 
  prev_saturation_liquid_key_ = Keys::getKey(domain_, "prev_saturation_liquid"); 

  // optional keys
  pressure_head_key_ = Keys::getKey(domain_, "pressure_head"); 

  // set up the base class 
  Flow_PK::Setup(S);

  // Our decision can be affected by the list of models
  Teuchos::RCP<Teuchos::ParameterList> physical_models =
      Teuchos::sublist(fp_list_, "physical models and assumptions");
  std::string mu_model = physical_models->get<std::string>("viscosity model", "constant viscosity");
  if (mu_model != "constant viscosity") {
    Errors::Message msg;
    msg << "Darcy PK supports only constant viscosity model.";
    Exceptions::amanzi_throw(msg);
  }
  // -- type of the flow (in matrix or on manifold)
  flow_on_manifold_ = physical_models->get<bool>("flow in fractures", false);
  flow_on_manifold_ &= (mesh_->manifold_dimension() != mesh_->space_dimension());

  // -- coupling with other PKs
  coupled_to_matrix_ = physical_models->get<std::string>("coupled matrix fracture flow", "") == "fracture";
  coupled_to_fracture_ = physical_models->get<std::string>("coupled matrix fracture flow", "") == "matrix";

  // Require primary field for this PK.
  Teuchos::RCP<Teuchos::ParameterList> list1 = Teuchos::sublist(fp_list_, "operators", true);
  Teuchos::RCP<Teuchos::ParameterList> list2 = Teuchos::sublist(list1, "diffusion operator", true);
  Teuchos::RCP<Teuchos::ParameterList> list3 = Teuchos::sublist(list2, "matrix", true);
  std::string name = list3->get<std::string>("discretization primary");

  if (!S->HasField(pressure_key_)) {
    std::vector<std::string> names;
    std::vector<AmanziMesh::Entity_kind> locations;
    std::vector<int> ndofs;

    names.push_back("cell");
    locations.push_back(AmanziMesh::CELL);
    ndofs.push_back(1);

    if (name != "fv: default" && name != "nlfv: default") {
      names.push_back("face");
      locations.push_back(AmanziMesh::FACE);
      ndofs.push_back(1);
    }

    S->RequireField(pressure_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", pressure_key_);
    auto eval = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S->SetFieldEvaluator(pressure_key_, eval);
  }

  // require additional fields for this PK
  if (!S->HasField(specific_storage_key_)) {
    S->RequireField(specific_storage_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField(specific_yield_key_)) {
    S->RequireField(specific_yield_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField(saturation_liquid_key_)) {
    S->RequireField(saturation_liquid_key_, saturation_liquid_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", saturation_liquid_key_);
    auto eval = Teuchos::rcp(new ConstantVariableFieldEvaluator(elist));
    S->SetFieldEvaluator(saturation_liquid_key_, eval);
  }

  if (!S->HasField(prev_saturation_liquid_key_)) {
    S->RequireField(prev_saturation_liquid_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField(darcy_flux_key_)) {
    if (flow_on_manifold_) {
      auto cvs = Operators::CreateNonManifoldCVS(mesh_);
      *S->RequireField(darcy_flux_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true) = *cvs;
    } else {
      S->RequireField(darcy_flux_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
        ->SetComponent("face", AmanziMesh::FACE, 1);
    }
  }

  {
    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", darcy_flux_key_);
    darcy_flux_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S->SetFieldEvaluator(darcy_flux_key_, darcy_flux_eval_);
  }

  // Require additional field evaluators for this PK.
  // -- porosity
  if (!S->HasField(porosity_key_)) {
    S->RequireField(porosity_key_, porosity_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(porosity_key_);
  }

  // -- viscosity
  if (!S->HasField("const_fluid_viscosity")) {
    S->RequireScalar("const_fluid_viscosity", passwd_);
  }

  // -- effective fracture permeability
  if (flow_on_manifold_) {
    if (!S->HasField(permeability_key_)) {
      S->RequireField(permeability_key_, permeability_key_)->SetMesh(mesh_)->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);

      Teuchos::RCP<Teuchos::ParameterList>
          fpm_list = Teuchos::sublist(fp_list_, "fracture permeability models", true);
      Teuchos::RCP<FracturePermModelPartition> fpm = CreateFracturePermModelPartition(mesh_, fpm_list);

      Teuchos::ParameterList elist;
      elist.set<std::string>("permeability key", permeability_key_)
           .set<std::string>("aperture key", Keys::getKey(domain_, "aperture"));
      Teuchos::RCP<FracturePermModelEvaluator> eval = Teuchos::rcp(new FracturePermModelEvaluator(elist, fpm));
      S->SetFieldEvaluator(permeability_key_, eval);
    }
  // -- matrix absolute permeability
  } else {
    if (!S->HasField(permeability_key_)) {
      S->RequireField(permeability_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, dim);
    }
  }

  // Local fields and evaluators.
  if (!S->HasField(hydraulic_head_key_)) {
    S->RequireField(hydraulic_head_key_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // full velocity vector
  if (!S->HasField(darcy_velocity_key_)) {
    S->RequireField(darcy_velocity_key_, darcy_velocity_key_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim);

    Teuchos::ParameterList elist;
    elist.set<std::string>("domain name", domain_);
    elist.set<std::string>("darcy velocity key", darcy_velocity_key_)
         .set<std::string>("darcy flux key", darcy_flux_key_);
    Teuchos::RCP<DarcyVelocityEvaluator> eval = Teuchos::rcp(new DarcyVelocityEvaluator(elist));
    S->SetFieldEvaluator(darcy_velocity_key_, eval);
  }

  // Require additional components for the existing fields
  Teuchos::ParameterList abs_perm = fp_list_->sublist("absolute permeability");
  coordinate_system_ = abs_perm.get<std::string>("coordinate system", "cartesian");
  int noff = abs_perm.get<int>("off-diagonal components", 0);
 
  if (noff > 0) {
    CompositeVectorSpace& cvs = *S->RequireField(permeability_key_, passwd_);
    cvs.SetOwned(false);
    cvs.AddComponent("offd", AmanziMesh::CELL, noff)->SetOwned(true);
  }

  // require optional fields
  if (fp_list_->isParameter("optional fields")) {
    std::vector<std::string> fields = fp_list_->get<Teuchos::Array<std::string> >("optional fields").toVector();
    for (auto it = fields.begin(); it != fields.end(); ++it) {
      Key optional_key = Keys::getKey(domain_, *it); 
      if (!S->HasField(optional_key)) {
        S->RequireField(optional_key, passwd_)->SetMesh(mesh_)->SetGhosted(true)
          ->SetComponent("cell", AmanziMesh::CELL, 1);
      }
    }
  }
}


/* ******************************************************************
* Extract information from parameter list and initialize data.
****************************************************************** */
void Darcy_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  // Initialize miscalleneous defaults.
  // -- times
  double t_ini = S->time(); 
  dt_next_ = dt_;
  dt_desirable_ = dt_;  // The minimum desirable time step from now on.
  dt_history_.clear();

  // -- others
  initialize_with_darcy_ = true;
  num_itrs_ = 0;

  // Create verbosity object to print out initialization statisticsr.,
  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = fp_list_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("DarcyPK-" + domain_, vlist)); 

  // Initilize various base class data.
  Flow_PK::Initialize(S);

  // Initialize local fields and evaluators. 
  InitializeFields_();
  UpdateLocalFields_(S);

  // Create solution and auxiliary data for time history.
  solution = S->GetFieldData(pressure_key_, passwd_);
  soln_->SetData(solution); 

  const Epetra_BlockMap& cmap = mesh_->cell_map(false);
  pdot_cells_prev = Teuchos::rcp(new Epetra_Vector(cmap));
  pdot_cells = Teuchos::rcp(new Epetra_Vector(cmap));
  
  std::string ti_method_name = ti_list_->get<std::string>("time integration method", "none");
  if (ti_method_name == "BDF1") {
    Teuchos::ParameterList& bdf1_list = ti_list_->sublist("BDF1");

    error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;  // usually 1e-4;

    // time step controller
    TimestepControllerFactory<Epetra_MultiVector> fac;
    ts_control_ = fac.Create(bdf1_list, pdot_cells, pdot_cells_prev);
  } else {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "WARNING: BDF1 time integration list is missing..." << std::endl;
  }

  // Initialize specific yield
  specific_yield_copy_ = Teuchos::null;
  UpdateSpecificYield_();

  // Initialize lambdas. It may be used by boundary conditions.
  CompositeVector& pressure = *S->GetFieldData(pressure_key_, passwd_);

  if (ti_list_->isSublist("pressure-lambda constraints") && solution->HasComponent("face")) {
    std::string method = ti_list_->sublist("pressure-lambda constraints").get<std::string>("method");
    if (method == "projection") {
      Epetra_MultiVector& p = *solution->ViewComponent("cell");
      Epetra_MultiVector& lambda = *solution->ViewComponent("face");

      DeriveFaceValuesFromCellValues(p, lambda);
    }
  }

  // Create and initialize boundary conditions and source terms.
  flux_units_ = 1.0;
  InitializeBCsSources_(*fp_list_);
  UpdateSourceBoundaryData(t_ini, t_ini, pressure);

  // Initialize diffusion operator and solver.
  // -- instead of scaling K, we scale the elemental mass matrices 
  double mu = *S->GetScalarData("const_fluid_viscosity");
  Teuchos::ParameterList& oplist = fp_list_->sublist("operators")
                                            .sublist("diffusion operator")
                                            .sublist("matrix");
  if (flow_on_manifold_)
      oplist.set<std::string>("nonlinear coefficient", "standard: cell");

  double factor = rho_ * rho_ / mu;
  if (coupled_to_fracture_) factor = rho_;

  Operators::PDE_DiffusionFactory opfactory;
  op_diff_ = opfactory.Create(oplist, mesh_, op_bc_, factor, gravity_);
  op_diff_->SetBCs(op_bc_, op_bc_);

  if (!flow_on_manifold_) {
    SetAbsolutePermeabilityTensor();
    Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K);
    op_diff_->Setup(Kptr, Teuchos::null, Teuchos::null);
  } else {
    S_->GetFieldEvaluator(permeability_key_)->HasFieldChanged(S_.ptr(), permeability_key_);
    auto Kptr = S_->GetFieldData(permeability_key_);
    op_diff_->Setup(Teuchos::null, Kptr, Teuchos::null);
  }

  op_diff_->ScaleMassMatrices(rho_ / mu);
  op_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op_ = op_diff_->global_operator();

  // -- accumulation operator.
  op_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op_));

  op_->SymbolicAssembleMatrix();
  op_->CreateCheckPoint();

  // -- generic linear solver.
  AMANZI_ASSERT(ti_list_->isParameter("linear solver"));
  solver_name_ = ti_list_->get<std::string>("linear solver");

  // -- preconditioner. There is no need to enhance it for Darcy
  AMANZI_ASSERT(ti_list_->isParameter("preconditioner"));
  std::string name = ti_list_->get<std::string>("preconditioner");
  Teuchos::ParameterList pc_list = preconditioner_list_->sublist(name);
  op_->InitializePreconditioner(pc_list);
  
  // Optional step: calculate hydrostatic solution consistent with BCs.
  // We have to do it only once per time period.
  bool init_darcy(false);
  if (ti_list_->isSublist("initialization") && initialize_with_darcy_) {
    initialize_with_darcy_ = false;
    bool wells_on = ti_list_->sublist("initialization").get<bool>("active wells", false);
    SolveFullySaturatedProblem(*solution, wells_on);
    init_darcy = true;
  }

  // Verbose output of initialization statistics.
  InitializeStatistics_(init_darcy);
}


/* ****************************************************************
* This completes initialization of common fields that were not 
* initialized by the state.
**************************************************************** */
void Darcy_PK::InitializeFields_()
{
  Teuchos::OSTab tab = vo_->getOSTab();

  InitializeField(S_.ptr(), passwd_, saturation_liquid_key_, 1.0);
  InitializeField(S_.ptr(), passwd_, prev_saturation_liquid_key_, 1.0);
}


/* ******************************************************************
* Print the header for new time period.
****************************************************************** */
void Darcy_PK::InitializeStatistics_(bool init_darcy)
{
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    double mu = *S_->GetScalarData("const_fluid_viscosity");
    std::string ti_name = ti_list_->get<std::string>("time integration method", "none");
    std::string ts_name = (ti_name == "BDF1") ? ti_list_->sublist(ti_name).get<std::string>("timestep controller type") 
                                              : "timestep controller is not defined";
    std::string pc_name = ti_list_->get<std::string>("preconditioner");

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "\nTI:\"" << ti_name.c_str() << "\""
               << " dt:" << ts_name
               << " LS:\"" << solver_name_.c_str() << "\""
               << " PC:\"" << pc_name.c_str() << "\"" << std::endl
               << "matrix: " << op_->PrintDiagnostics() << std::endl;
    *vo_->os() << "constant viscosity model, mu=" << mu << std::endl;

    if (init_darcy) {
      *vo_->os() << "initial pressure guess: \"from saturated solver\"\n" << std::endl;
    } else {
      *vo_->os() << "initial pressure guess: \"from State\"\n" << std::endl;
    }

    int missed_tmp = missed_bc_faces_;
    int dirichlet_tmp = dirichlet_bc_faces_;
#ifdef HAVE_MPI
    mesh_->get_comm()->SumAll(&missed_tmp, &missed_bc_faces_, 1);
    mesh_->get_comm()->SumAll(&dirichlet_tmp, &dirichlet_bc_faces_, 1);
#endif

    *vo_->os() << "pressure BC assigned to " << dirichlet_bc_faces_ << " faces" << std::endl;
    *vo_->os() << "default (no-flow) BC assigned to " << missed_bc_faces_ << " faces" << std::endl << std::endl;

    VV_PrintHeadExtrema(*solution);
    VV_PrintSourceExtrema();

    *vo_->os() << vo_->color("green") << "Initialization of PK is complete, T=" 
               << S_->time() << " dT=" << dt_ << vo_->reset() << std::endl << std::endl;
  }
}


/* ******************************************************************* 
* Performs one time step from t_old to t_new. The boundary conditions
* are calculated only once, during the initialization step.  
******************************************************************* */
bool Darcy_PK::AdvanceStep(double t_old, double t_new, bool reinit) 
{
  dt_ = t_new - t_old;
  double dt_MPC(dt_);

  // refresh data
  UpdateSourceBoundaryData(t_old, t_new, *solution);

  // calculate and assemble elemental stiffness matrices
  double factor = 1.0 / g_;
  const CompositeVector& ss = *S_->GetFieldData(specific_storage_key_);
  CompositeVector ss_g(ss); 
  ss_g.Update(0.0, ss, factor);

  factor = 1.0 / (g_ * dt_);
  CompositeVector sy_g(*specific_yield_copy_); 
  sy_g.Scale(factor);

  op_->RestoreCheckPoint();
  op_acc_->AddAccumulationDelta(*solution, ss_g, ss_g, dt_, "cell");
  op_acc_->AddAccumulationDeltaNoVolume(*solution, sy_g, "cell");

  // Peaceman model
  if (S_->HasField("well_index")) {
    const CompositeVector& wi = *S_->GetFieldData("well_index");
    op_acc_->AddAccumulationTerm(wi, "cell");
  }

  op_diff_->ApplyBCs(true, true, true);

  CompositeVector& rhs = *op_->rhs();
  AddSourceTerms(rhs);

  op_->AssembleMatrix();
  op_->UpdatePreconditioner();

  // save pressure at time t^n.
  std::string dt_control = ti_list_->sublist("BDF1").get<std::string>("timestep controller type");
  Teuchos::RCP<Epetra_MultiVector> p_old;
  if (dt_control == "adaptive") {
    p_old = Teuchos::rcp(new Epetra_MultiVector(*solution->ViewComponent("cell")));
  }

  // create linear solver and calculate new pressure
  AmanziSolvers::LinearOperatorFactory<Operators::Operator, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::Operator, CompositeVector, CompositeVectorSpace> >
     solver = factory.Create(solver_name_, *linear_operator_list_, op_);

  solver->add_criteria(AmanziSolvers::LIN_SOLVER_MAKE_ONE_ITERATION);
  solver->ApplyInverse(rhs, *solution);

  // statistics
  num_itrs_++;

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    double pnorm;
    solution->Norm2(&pnorm);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "pressure solver (" << solver->name()
               << "): ||p,lambda||=" << pnorm 
               << "  itrs=" << solver->num_itrs() << std::endl;
    VV_PrintHeadExtrema(*solution);
  }

  // calculate time derivative and 2nd-order solution approximation
  if (dt_control == "adaptive") {
    Epetra_MultiVector& p_new = *solution->ViewComponent("cell");  // pressure at t^{n+1}

    for (int c = 0; c < ncells_owned; c++) {
      (*pdot_cells)[c] = (p_new[0][c] - (*p_old)[0][c]) / dt_; 
      p_new[0][c] = (*p_old)[0][c] + ((*pdot_cells_prev)[c] + (*pdot_cells)[c]) * dt_ / 2;
    }
  }

  // estimate time multiplier
  dt_desirable_ = ts_control_->get_timestep(dt_MPC, 1);

  // Darcy_PK always takes the suggested time step and cannot fail
  dt_tuple times(t_new, dt_MPC);
  dt_history_.push_back(times);

  return false;
}


/* ******************************************************************
* Transfer data from the external flow state FS_MPC. MPC may request
* to populate the original state FS. 
****************************************************************** */
void Darcy_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  // calculate Darcy mass flux
  Teuchos::RCP<CompositeVector> flux = S->GetFieldData(darcy_flux_key_, passwd_);

  if (coupled_to_matrix_ || flow_on_manifold_) {
    op_diff_->UpdateFluxNonManifold(solution.ptr(), flux.ptr());
    flux->Scale(1.0 / rho_);
    FractureConservationLaw_();
  } else {
    op_diff_->UpdateFlux(solution.ptr(), flux.ptr());
    flux->Scale(1.0 / rho_);
  }

  // update time derivative
  *pdot_cells_prev = *pdot_cells;
}


/* ******************************************************************
* Add area/length factor to specific yield. 
****************************************************************** */
void Darcy_PK::UpdateSpecificYield_()
{
  specific_yield_copy_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData(specific_yield_key_), true));

  // do we have non-zero specific yield? 
  double tmp;
  specific_yield_copy_->Norm2(&tmp);
  if (tmp == 0.0) return;

  // populate ghost cells
  specific_yield_copy_->ScatterMasterToGhosted();
  const Epetra_MultiVector& specific_yield = *specific_yield_copy_->ViewComponent("cell", true);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  int negative_yield = 0;
  for (int c = 0; c < ncells_owned; c++) {
    if (specific_yield[0][c] > 0.0) {
      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

      double area = 0.0;
      int nfaces = faces.size();
      for (int n = 0; n < nfaces; n++) {
        int f = faces[n];
        int c2 = WhetStone::cell_get_face_adj_cell(*mesh_, c, f);

        if (c2 >= 0) {
          if (specific_yield[0][c2] <= 0.0)  // cell in the fully saturated layer
            area -= (mesh_->face_normal(f))[dim - 1] * dirs[n];
        }
      }
      specific_yield[0][c] *= area;
      if (area <= 0.0) negative_yield++;
    }
  }

#ifdef HAVE_MPI
  int negative_yield_tmp = negative_yield;
  mesh_->get_comm()->MaxAll(&negative_yield_tmp, &negative_yield, 1);
#endif
  if (negative_yield > 0) {
    Errors::Message msg;
    msg << "Flow PK: configuration of the yield region leads to negative yield interfaces.";
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
* This is strange.
****************************************************************** */
void Darcy_PK::CalculateDiagnostics(const Teuchos::RCP<State>& S) {
  UpdateLocalFields_(S.ptr());
}


/* ******************************************************************
* TBW
****************************************************************** */
void Darcy_PK::FractureConservationLaw_()
{
  if (!coupled_to_matrix_) return;

  const auto& fracture_flux = *S_->GetFieldData("fracture-darcy_flux")->ViewComponent("face", true);
  const auto& matrix_flux = *S_->GetFieldData("darcy_flux")->ViewComponent("face", true);

  const auto& fracture_map = S_->GetFieldData("fracture-darcy_flux")->Map().Map("face", true);
  const auto& matrix_map = S_->GetFieldData("darcy_flux")->Map().Map("face", true);

  auto mesh_matrix = S_->GetMesh("domain");

  std::vector<int> dirs;
  double err(0.0), flux_max(0.0);
  AmanziMesh::Entity_ID_List faces, cells;

  for (int c = 0; c < ncells_owned; c++) {
    double flux_sum(0.0);
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    for (int i = 0; i < faces.size(); i++) {
      int f = faces[i];
      int g = fracture_map->FirstPointInElement(f);
      int ndofs = fracture_map->ElementSize(f);
      if (ndofs > 1) g += Operators::UniqueIndexFaceToCells(*mesh_, f, c);

      flux_sum += fracture_flux[0][g] * dirs[i];
    }

    // sum into fluxes from matrix
    auto f = mesh_->entity_get_parent(AmanziMesh::CELL, c);
    mesh_matrix->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int pos = Operators::UniqueIndexFaceToCells(*mesh_matrix, f, cells[0]);

    for (int j = 0; j != cells.size(); ++j) {
      mesh_matrix->cell_get_faces_and_dirs(cells[j], &faces, &dirs);

      for (int i = 0; i < faces.size(); i++) {
        if (f == faces[i]) {
          int g = matrix_map->FirstPointInElement(f);            
          double fln = matrix_flux[0][g + (pos + j) % 2] * dirs[i];
          flux_sum -= fln;
          flux_max = std::max<double>(flux_max, std::fabs(fln));
            
          break;
        }
      }
    }

    err = std::max<double>(err, fabs(flux_sum));
  }

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "maximum error in conservation law:" << err << " flux_max=" << flux_max << std::endl;
  }
}

}  // namespace Flow
}  // namespace Amanzi

