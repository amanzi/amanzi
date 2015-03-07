/*
  This is the flow component of the Amanzi code.

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

#include "errors.hh"
#include "exceptions.hh"
#include "LinearOperatorFactory.hh"
#include "mfd3d_diffusion.hh"
#include "OperatorDiffusionFactory.hh"
#include "tensor.hh"

#include "Darcy_PK.hh"
#include "FlowDefs.hh"
#include "Flow_SourceFactory.hh"

#include "DarcyVelocityEvaluator.hh"
#include "primary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Simplest possible constructor: extracts lists and requires fields.
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
  preconditioner_list_ = Teuchos::sublist(glist, "Preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "Solvers", true);
  ti_list_ = Teuchos::sublist(dp_list_, "time integrator", true);
}


/* ******************************************************************
* Define structure of this PK.
****************************************************************** */
void Darcy_PK::Setup()
{
  mesh_ = S_->GetMesh();
  dim = mesh_->space_dimension();

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
  if (name != "fv: default") {
    names.push_back("face");
    locations.push_back(AmanziMesh::FACE);
    ndofs.push_back(1);
  }

  if (!S_->HasField("pressure")) {
    S_->RequireField("pressure", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);
  }

  // require additional fields for this PK
  if (!S_->HasField("fluid_density")) {
    S_->RequireScalar("fluid_density", passwd_);
  }
  if (!S_->HasField("fluid_viscosity")) {
    S_->RequireScalar("fluid_viscosity", passwd_);
  }
  if (!S_->HasField("gravity")) {
    S_->RequireConstantVector("gravity", passwd_, dim);  // state resets ownership.
  } 

  if (!S_->HasField("permeability")) {
    S_->RequireField("permeability", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim);
  }

  if (!S_->HasField("specific_storage")) {
    S_->RequireField("specific_storage", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S_->HasField("specific_yield")) {
    S_->RequireField("specific_yield", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S_->HasField("saturation_liquid")) {
    S_->RequireField("saturation_liquid", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S_->HasField("prev_saturation_liquid")) {
    S_->RequireField("prev_saturation_liquid", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S_->HasField("darcy_flux")) {
    S_->RequireField("darcy_flux", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);

    Teuchos::ParameterList elist;
    elist.set<std::string>("evaluator name", "darcy_flux");
    darcy_flux_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
    S_->SetFieldEvaluator("darcy_flux", darcy_flux_eval_);
  }

  // Require additional field evaluators for this PK.
  // porosity
  if (!S_->HasField("porosity")) {
    S_->RequireField("porosity", "porosity")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireFieldEvaluator("porosity");
  }

  // Local fields and evaluators.
  if (!S_->HasField("hydraulic_head")) {
    S_->RequireField("hydraulic_head", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // full velocity vector
  if (!S_->HasField("darcy_velocity")) {
    S_->RequireField("darcy_velocity", "darcy_velocity")->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim);

    Teuchos::ParameterList elist;
    Teuchos::RCP<DarcyVelocityEvaluator> eval = Teuchos::rcp(new DarcyVelocityEvaluator(elist));
    S_->SetFieldEvaluator("darcy_velocity", eval);
  }
}


/* ******************************************************************
* Clean memory.
****************************************************************** */
Darcy_PK::~Darcy_PK()
{
  if (bc_pressure != NULL) delete bc_pressure;
  if (bc_head != NULL) delete bc_head;
  if (bc_flux != NULL) delete bc_flux;
  if (bc_seepage != NULL) delete bc_seepage;

  if (src_sink != NULL) delete src_sink;
  if (vo_ != NULL) delete vo_;
}


/* ******************************************************************
* Extract information from Diffusion problem parameter list.
****************************************************************** */
void Darcy_PK::Initialize()
{
  // Initialize defaults
  bc_seepage = NULL; 
  src_sink = NULL;
  src_sink_distribution = 0;

  initialize_with_darcy_ = true;
  num_itrs_ = 0;

  // Time control specific to this PK.
  ResetPKtimes(0.0, FLOW_INITIAL_DT);
  dT_desirable_ = dT;

  // create verbosity object
  Teuchos::ParameterList vlist;
  vlist.sublist("VerboseObject") = dp_list_->sublist("VerboseObject");
  vo_ = new VerboseObject("FlowPK::Darcy", vlist); 

  // Initilize various common data depending on mesh and state.
  Flow_PK::Initialize();

  // Create local evaluators. Initialize local fields.
  InitializeFields_();
  UpdateLocalFields_();

  // Initialize BCs and source terms.
  InitializeBCsSources_(*dp_list_);
  op_bc_ = Teuchos::rcp(new Operators::BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model, bc_value, bc_mixed));

  // Create solution and auxiliary data for time history.
  // solution = Teuchos::rcp(new CompositeVector(*(S_->GetFieldData("pressure"))));
  solution = S_->GetFieldData("pressure", passwd_);

  const Epetra_BlockMap& cmap = mesh_->cell_map(false);
  pdot_cells_prev = Teuchos::rcp(new Epetra_Vector(cmap));
  pdot_cells = Teuchos::rcp(new Epetra_Vector(cmap));
  
  // Initialize times.
  double time = S_->time();
  if (time >= 0.0) T_physics = time;

  // Initialize boundary condtions. 
  ProcessShiftWaterTableList(*dp_list_);

  time = T_physics;
  bc_pressure->Compute(time);
  bc_flux->Compute(time);
  bc_seepage->Compute(time);
  if (shift_water_table_.getRawPtr() == NULL) {
    bc_head->Compute(time);
  } else {
    bc_head->ComputeShift(time, shift_water_table_->Values());
  }

  const CompositeVector& pressure = *S_->GetFieldData("pressure");
  ComputeBCs(pressure);

  // Allocate memory for other fundamental structures
  K.resize(ncells_owned);

  if (src_sink_distribution & CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
    Kxy = Teuchos::rcp(new Epetra_Vector(mesh_->cell_map(true)));
  }
}


/* ******************************************************************
* Initialization of auxiliary variables (lambda and two saturations).
* WARNING: Flow_PK may use complex initialization of the remaining 
* state variables.
****************************************************************** */
void Darcy_PK::InitializeAuxiliaryData()
{
  // pressures (lambda is not important when solver is very accurate)
  CompositeVector& cv = *S_->GetFieldData("pressure", passwd_);
  const Epetra_MultiVector& pressure = *(cv.ViewComponent("cell"));
  Epetra_MultiVector& lambda = *(cv.ViewComponent("face"));

  DeriveFaceValuesFromCellValues(pressure, lambda);

  // saturations
  if (!S_->GetField("saturation_liquid", passwd_)->initialized()) {
    S_->GetFieldData("saturation_liquid", passwd_)->PutScalar(1.0);
    S_->GetField("saturation_liquid", passwd_)->set_initialized();
  }
  if (!S_->GetField("prev_saturation_liquid", passwd_)->initialized()) {
    S_->GetFieldData("prev_saturation_liquid", passwd_)->PutScalar(1.0);
    S_->GetField("prev_saturation_liquid", passwd_)->set_initialized();
  }
}


/* ******************************************************************
* Specific initialization of a steady state time integration phase.
* WARNING: now it is equivalent to transient phase.
****************************************************************** */
void Darcy_PK::InitTimeInterval()
{
  // times
  double T0 = ti_list_->get<double>("start interval time", 0.0);
  double dT0 = ti_list_->get<double>("initial time step", 1.0);
  ResetPKtimes(T0, dT0);

  dT = dT0;
  dTnext = dT0;
  dT_desirable_ = dT0;  // The minimum desirable time step from now on.
  dT_history_.clear();

  std::string ti_method_name = ti_list_->get<std::string>("time integration method", "none");
  ASSERT(ti_method_name == "BDF1");
  Teuchos::ParameterList& bdf1_list = ti_list_->sublist("BDF1");

  std::string dT_method_name = bdf1_list.get<std::string>("timestep controller type");
  Teuchos::ParameterList dtlist;
  if (dT_method_name == "standard") {
    dtlist = bdf1_list.sublist("timestep controller standard parameters");
    dTfactor_ = dtlist.get<double>("time step increase factor");
  } else if (dT_method_name == "fixed") {
    dtlist = bdf1_list.sublist("timestep controller fixed parameters");
    dTfactor_ = dtlist.get<double>("time step increase factor");
  } else if (dT_method_name == "adaptive") {
    dtlist = bdf1_list.sublist("timestep controller adaptive parameters");
  }
  dTmax_ = dtlist.get<double>("max time step", Flow::FLOW_MAXIMUM_DT);

  // initialize error control
  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;  // usually 1e-4;

  // initialize specific yield
  specific_yield_copy_ = Teuchos::null;
  UpdateSpecificYield_();

  // initialize diffusion operator
  SetAbsolutePermeabilityTensor();

  Teuchos::ParameterList& oplist = dp_list_->sublist("operators")
                                            .sublist("diffusion operator")
                                            .sublist("matrix");
  Operators::OperatorDiffusionFactory opfactory;
  op_diff_ = opfactory.Create(mesh_, op_bc_, oplist, gravity_, 0);  // The last 0 means no upwind
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K);
  op_diff_->SetBCs(op_bc_);
  op_diff_->Setup(Kptr, Teuchos::null, Teuchos::null, rho_, mu_);
  op_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op_ = op_diff_->global_operator();

  // initialize accumulation operator
  op_acc_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, op_));

  op_->SymbolicAssembleMatrix();
  op_->CreateCheckPoint();

  // preconditioner and optional linear solver
  ASSERT(ti_list_->isParameter("preconditioner"));
  preconditioner_name_ = ti_list_->get<std::string>("preconditioner");
  ASSERT(preconditioner_list_->isSublist(preconditioner_name_));
  
  ASSERT(ti_list_->isParameter("linear solver"));
  solver_name_ = ti_list_->get<std::string>("linear solver");

  // initialize well modeling
  if (src_sink != NULL) {
    double T1 = T0 + dT0;
    if (src_sink_distribution & CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      CalculatePermeabilityFactorInWell();
      src_sink->ComputeDistribute(T0, T1, Kxy->Values()); 
    } else {
      src_sink->ComputeDistribute(T0, T1, NULL);
    }
  }
  
  // Optional step: calculate hydrostatic solution consistent with BCs.
  // We have to do it only once per time period.
  if (ti_list_->isSublist("initialization") && initialize_with_darcy_) {
    initialize_with_darcy_ = false;
    Epetra_MultiVector& p = *solution->ViewComponent("cell");
    Epetra_MultiVector& lambda = *solution->ViewComponent("face");
    DeriveFaceValuesFromCellValues(p, lambda);

    SolveFullySaturatedProblem(T0, *solution);
    pressure_eval_->SetFieldAsChanged(S_.ptr());
  }

  // print initialization head for this time period
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << std::endl 
        << vo_->color("green") << "New TI period: " << ti_method_name.c_str() 
        << vo_->reset() << std::endl;
    *vo_->os()<< "dT controller:" << dT_method_name << " Src:" << src_sink_distribution
              << " LS:\"" << solver_name_.c_str() << "\""
              << " PC:\"" << preconditioner_name_.c_str() << "\"" << std::endl;

    if (initialize_with_darcy_) {
      *vo_->os() << "initial pressure guess: \"saturated solution\"\n" << std::endl;
    } else {
      *vo_->os() << "initial pressure guess: \"from State\"\n" << std::endl;
    }

    VV_PrintHeadExtrema(*solution);
  }
}


/* ******************************************************************* 
* Performs one time step of size dT. The boundary conditions are 
* calculated only once, during the initialization step.  
******************************************************************* */
bool Darcy_PK::Advance(double dT_MPC, double& dT_actual) 
{
  dT = dT_MPC;
  double T1 = S_->time();
  if (T1 >= 0.0) T_physics = T1;

  // update boundary conditions and source terms
  T1 = T_physics;
  bc_pressure->Compute(T1);
  bc_flux->Compute(T1);
  bc_seepage->Compute(T1);
  if (shift_water_table_.getRawPtr() == NULL)
    bc_head->Compute(T1);
  else
    bc_head->ComputeShift(T1, shift_water_table_->Values());

  if (src_sink != NULL) {
    double T0 = T1 - dT_MPC; 
    if (src_sink_distribution & CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      src_sink->ComputeDistribute(T0, T1, Kxy->Values()); 
    } else {
      src_sink->ComputeDistribute(T0, T1, NULL);
    }
  }

  ComputeBCs(*solution);

  // calculate and assemble elemental stifness matrices
  double factor = 1.0 / g_;
  const CompositeVector& ss = *S_->GetFieldData("specific_storage");
  CompositeVector ss_g(ss); 
  ss_g.Update(0.0, ss, factor);

  factor = 1.0 / (g_ * dT);
  CompositeVector sy_g(*specific_yield_copy_); 
  sy_g.Scale(factor);

  op_->RestoreCheckPoint();
  op_acc_->AddAccumulationTerm(*solution, ss_g, dT, "cell");
  op_acc_->AddAccumulationTerm(*solution, sy_g, "cell");

  op_diff_->ApplyBCs(true);
  op_->AssembleMatrix();
  op_->InitPreconditioner(preconditioner_name_, *preconditioner_list_);

  CompositeVector& rhs = *op_->rhs();
  if (src_sink != NULL) AddSourceTerms(rhs);

  // create linear solver
  AmanziSolvers::LinearOperatorFactory<Operators::Operator, CompositeVector, CompositeVectorSpace> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::Operator, CompositeVector, CompositeVectorSpace> >
     solver = factory.Create(solver_name_, *linear_operator_list_, op_);

  solver->add_criteria(AmanziSolvers::LIN_SOLVER_MAKE_ONE_ITERATION);
  solver->ApplyInverse(rhs, *solution);

  bool fail = false;

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
  std::string dT_method_name = ti_list_->sublist("BDF1").get<std::string>("timestep controller type");

  if (dT_method_name == "adaptive") {
    const Epetra_MultiVector& p = *S_->GetFieldData("pressure")->ViewComponent("cell");  // pressure at t^n
    Epetra_MultiVector& p_cell = *solution->ViewComponent("cell");  // pressure at t^{n+1}

    for (int c = 0; c < ncells_owned; c++) {
      (*pdot_cells)[c] = (p_cell[0][c] - p[0][c]) / dT; 
      p_cell[0][c] = p[0][c] + ((*pdot_cells_prev)[c] + (*pdot_cells)[c]) * dT / 2;
    }
  }

  // estimate time multiplier
  if (dT_method_name == "adaptive") {
    double err, dTfactor;
    err = ErrorEstimate_(&dTfactor);
    if (err > 0.0) throw 1000;  // fix (lipnikov@lan.gov)
    dT_desirable_ = std::min(dT_MPC * dTfactor_, dTmax_);
  } else {
    dT_desirable_ = std::min(dT_desirable_ * dTfactor_, dTmax_);
  }

  // Darcy_PK always takes suggested time step
  dT_actual = dT_MPC;

  dt_tuple times(T1, dT_MPC);
  dT_history_.push_back(times);

  return fail;
}


/* ******************************************************************
* Transfer data from the external flow state FS_MPC. MPC may request
* to populate the original state FS. 
****************************************************************** */
void Darcy_PK::CommitStep(double dt, const Teuchos::Ptr<State>& S)
{
  CompositeVector& p = *S->GetFieldData("pressure", passwd_);
  p = *solution;

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
void Darcy_PK::CalculateDiagnostics(const Teuchos::Ptr<State>& S) {
  UpdateLocalFields_();
}

}  // namespace Flow
}  // namespace Amanzi

