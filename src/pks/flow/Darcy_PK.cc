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

#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "boost/algorithm/string.hpp"

#include "DiffusionFactory.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "LinearOperatorFactory.hh"
#include "mfd3d_diffusion.hh"
#include "TimestepControllerFactory.hh"
#include "Tensor.hh"

#include "Darcy_PK.hh"
#include "FlowDefs.hh"

#include "DarcyVelocityEvaluator.hh"
#include "primary_variable_field_evaluator.hh"

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
  const char* result = pk_name.data();

  boost::iterator_range<std::string::iterator> res = boost::algorithm::find_last(pk_name,"->"); 
  if (res.end() - pk_name.end() != 0) boost::algorithm::erase_head(pk_name,  res.end() - pk_name.begin());

  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  Teuchos::RCP<Teuchos::ParameterList> flow_list = Teuchos::sublist(pk_list, pk_name, true);
  dp_list_ = Teuchos::sublist(flow_list, "Darcy problem", true);

  // We also need iscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(dp_list_, "time integrator", true);
}


/* ******************************************************************
* Old constructor for unit tests.
****************************************************************** */
Darcy_PK::Darcy_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                   const std::string& pk_list_name,
                   Teuchos::RCP<State> S) :
    Flow_PK()
{
  S_ = S;

  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  Teuchos::RCP<Teuchos::ParameterList> flow_list = Teuchos::sublist(pk_list, pk_list_name, true);
  dp_list_ = Teuchos::sublist(flow_list, "Darcy problem", true);

  // We also need iscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "solvers", true);
  ti_list_ = Teuchos::sublist(dp_list_, "time integrator", true);
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
  mesh_ = S->GetMesh();
  dim = mesh_->space_dimension();

  Flow_PK::Setup(S);

  // Our decision can be affected by the list of models
  Teuchos::RCP<Teuchos::ParameterList> physical_models =
      Teuchos::sublist(dp_list_, "physical models and assumptions");
  std::string mu_model = physical_models->get<std::string>("viscosity model", "constant viscosity");
  if (mu_model != "constant viscosity") {
    Errors::Message msg;
    msg << "Darcy PK supports only constant viscosity model.";
    Exceptions::amanzi_throw(msg);
  }

  // Require primary field for this PK.
  std::vector<std::string> names;
  std::vector<AmanziMesh::Entity_kind> locations;
  std::vector<int> ndofs;

  Teuchos::RCP<Teuchos::ParameterList> list1 = Teuchos::sublist(dp_list_, "operators", true);
  Teuchos::RCP<Teuchos::ParameterList> list2 = Teuchos::sublist(list1, "diffusion operator", true);
  Teuchos::RCP<Teuchos::ParameterList> list3 = Teuchos::sublist(list2, "matrix", true);
  std::string name = list3->get<std::string>("discretization primary");

  names.push_back("cell");
  locations.push_back(AmanziMesh::CELL);
  ndofs.push_back(1);
  if (name != "fv: default" && name != "nlfv: default") {
    names.push_back("face");
    locations.push_back(AmanziMesh::FACE);
    ndofs.push_back(1);
  }

  if (!S->HasField("pressure")) {
    S->RequireField("pressure", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);
  }

  // require additional fields for this PK
  if (!S->HasField("specific_storage")) {
    S->RequireField("specific_storage", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S->HasField("specific_yield")) {
    S->RequireField("specific_yield", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField("saturation_liquid")) {
    S->RequireField("saturation_liquid", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S->HasField("prev_saturation_liquid")) {
    S->RequireField("prev_saturation_liquid", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField("darcy_flux")) {
    S->RequireField("darcy_flux", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", "darcy_flux");
    darcy_flux_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S->SetFieldEvaluator("darcy_flux", darcy_flux_eval_);
  }

  // Require additional field evaluators for this PK.
  // -- porosity
  if (!S->HasField("porosity")) {
    S->RequireField("porosity", "porosity")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator("porosity");
  }

  // -- viscosity
  if (!S->HasField("fluid_viscosity")) {
    S->RequireScalar("fluid_viscosity", passwd_);
  }

  // Local fields and evaluators.
  if (!S->HasField("hydraulic_head")) {
    S->RequireField("hydraulic_head", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // full velocity vector
  if (!S->HasField("darcy_velocity")) {
    S->RequireField("darcy_velocity", "darcy_velocity")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim);

    Teuchos::ParameterList elist;
    Teuchos::RCP<DarcyVelocityEvaluator> eval = Teuchos::rcp(new DarcyVelocityEvaluator(elist));
    S->SetFieldEvaluator("darcy_velocity", eval);
  }

  // Require additional components for the existing fields
  Teuchos::ParameterList abs_perm = dp_list_->sublist("absolute permeability");
  coordinate_system_ = abs_perm.get<std::string>("coordinate system", "cartesian");
  int noff = abs_perm.get<int>("off-diagonal components", 0);
 
  if (noff > 0) {
    CompositeVectorSpace& cvs = *S->RequireField("permeability", passwd_);
    cvs.SetOwned(false);
    cvs.AddComponent("offd", AmanziMesh::CELL, noff)->SetOwned(true);
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
  vlist.sublist("verbose object") = dp_list_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("FlowPK::Darcy", vlist)); 

  // Initilize various base class data.
  Flow_PK::Initialize(S);

  // Initialize local fields and evaluators. 
  InitializeFields_();
  UpdateLocalFields_(S);

  // Create solution and auxiliary data for time history.
  solution = S->GetFieldData("pressure", passwd_);

  const Epetra_BlockMap& cmap = mesh_->cell_map(false);
  pdot_cells_prev = Teuchos::rcp(new Epetra_Vector(cmap));
  pdot_cells = Teuchos::rcp(new Epetra_Vector(cmap));
  
  std::string ti_method_name = ti_list_->get<std::string>("time integration method", "none");
  ASSERT(ti_method_name == "BDF1");
  Teuchos::ParameterList& bdf1_list = ti_list_->sublist("BDF1");

  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;  // usually 1e-4;

  // time step controller
  TimestepControllerFactory<Epetra_MultiVector> fac;
  ts_control_ = fac.Create(bdf1_list, pdot_cells, pdot_cells_prev);

  // Initialize specific yield
  specific_yield_copy_ = Teuchos::null;
  UpdateSpecificYield_();

  // Initialize lambdas. It may be used by boundary conditions.
  CompositeVector& pressure = *S->GetFieldData("pressure", passwd_);

  if (pressure.HasComponent("face")) {
    Epetra_MultiVector& p = *solution->ViewComponent("cell");
    Epetra_MultiVector& lambda = *solution->ViewComponent("face");

    DeriveFaceValuesFromCellValues(p, lambda);
  }

  // Create and initialize boundary conditions and source terms.
  flux_units_ = 1.0;
  InitializeBCsSources_(*dp_list_);
  UpdateSourceBoundaryData(t_ini, t_ini, pressure);

  // Initialize diffusion operator and solver.
  // -- instead of scaling K, we scale the elemental mass matrices 
  double mu = *S->GetScalarData("fluid_viscosity");
  SetAbsolutePermeabilityTensor();

  Teuchos::ParameterList& oplist = dp_list_->sublist("operators")
                                            .sublist("diffusion operator")
                                            .sublist("matrix");
  Operators::DiffusionFactory opfactory;
  op_diff_ = opfactory.Create(oplist, mesh_, op_bc_, rho_ * rho_ / mu, gravity_);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K);
  op_diff_->SetBCs(op_bc_, op_bc_);
  op_diff_->Setup(Kptr, Teuchos::null, Teuchos::null);
  op_diff_->ScaleMassMatrices(rho_ / mu);
  op_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op_ = op_diff_->global_operator();

  // -- accumulation operator.
  op_acc_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, op_));

  op_->SymbolicAssembleMatrix();
  op_->CreateCheckPoint();

  // -- generic linear solver.
  ASSERT(ti_list_->isParameter("linear solver"));
  solver_name_ = ti_list_->get<std::string>("linear solver");

  // -- preconditioner. There is no need to enhance it for Darcy
  ASSERT(ti_list_->isParameter("preconditioner"));
  preconditioner_name_ = ti_list_->get<std::string>("preconditioner");
  ASSERT(preconditioner_list_->isSublist(preconditioner_name_));
  
  // Optional step: calculate hydrostatic solution consistent with BCs.
  // We have to do it only once per time period.
  bool init_darcy(false);
  if (ti_list_->isSublist("initialization") && initialize_with_darcy_) {
    initialize_with_darcy_ = false;
    SolveFullySaturatedProblem(*solution);
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

  // set popular default values for missed fields.
  if (S_->GetField("saturation_liquid")->owner() == passwd_) {
    if (S_->HasField("saturation_liquid")) {
      if (!S_->GetField("saturation_liquid", passwd_)->initialized()) {
        S_->GetFieldData("saturation_liquid", passwd_)->PutScalar(1.0);
        S_->GetField("saturation_liquid", passwd_)->set_initialized();

        if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
            *vo_->os() << "initilized saturation_liquid to default value 1.0" << std::endl;  
      }
    }
  }

  if (S_->HasField("prev_saturation_liquid")) {
    if (!S_->GetField("prev_saturation_liquid", passwd_)->initialized()) {
      S_->GetFieldData("prev_saturation_liquid", passwd_)->PutScalar(1.0);
      S_->GetField("prev_saturation_liquid", passwd_)->set_initialized();

      if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM)
          *vo_->os() << "initilized prev_saturation_liquid to default value 1.0" << std::endl;  
    }
  }
}


/* ******************************************************************
* Print the header for new time period.
****************************************************************** */
void Darcy_PK::InitializeStatistics_(bool init_darcy)
{
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    double mu = *S_->GetScalarData("fluid_viscosity");
    std::string ti_method_name = ti_list_->get<std::string>("time integration method");
    std::string dt_method_name = ti_list_->sublist("BDF1").get<std::string>("timestep controller type");

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "\nTI:\"" << ti_method_name.c_str() << "\""
               << " dt:" << dt_method_name
               << " LS:\"" << solver_name_.c_str() << "\""
               << " PC:\"" << preconditioner_name_.c_str() << "\"" << std::endl
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

    *vo_->os() << vo_->color("green") << "Initalization of PK is complete, T=" 
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
  const CompositeVector& ss = *S_->GetFieldData("specific_storage");
  CompositeVector ss_g(ss); 
  ss_g.Update(0.0, ss, factor);

  factor = 1.0 / (g_ * dt_);
  CompositeVector sy_g(*specific_yield_copy_); 
  sy_g.Scale(factor);

  op_->RestoreCheckPoint();
  op_acc_->AddAccumulationTerm(*solution, ss_g, dt_, "cell");
  op_acc_->AddAccumulationTerm(*solution, sy_g, "cell");

  op_diff_->ApplyBCs(true, true);
  op_->AssembleMatrix();
  op_->InitPreconditioner(preconditioner_name_, *preconditioner_list_);

  CompositeVector& rhs = *op_->rhs();
  AddSourceTerms(rhs);

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
               << "): ||p,lambda||=" << pnorm << std::endl;
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
  // calculate darcy mass flux
  CompositeVector& darcy_flux = *S->GetFieldData("darcy_flux", passwd_);
  op_diff_->UpdateFlux(*solution, darcy_flux);

  Epetra_MultiVector& flux = *darcy_flux.ViewComponent("face", true);
  for (int f = 0; f < nfaces_owned; f++) flux[0][f] /= rho_;

  // update time derivative
  *pdot_cells_prev = *pdot_cells;
}


/* ******************************************************************
* Add area/length factor to specific yield. 
****************************************************************** */
void Darcy_PK::UpdateSpecificYield_()
{
  specific_yield_copy_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("specific_yield"), true));

  // do we have non-zero specific yield? 
  double tmp;
  specific_yield_copy_->Norm2(&tmp);
  if (tmp == 0.0) return;

  // populate ghost cells
  specific_yield_copy_->ScatterMasterToGhosted();
  const Epetra_MultiVector& specific_yield = *specific_yield_copy_->ViewComponent("cell", true);

  WhetStone::MFD3D_Diffusion mfd3d(mesh_);
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
        int c2 = mfd3d.cell_get_face_adj_cell(c, f);

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

}  // namespace Flow
}  // namespace Amanzi

