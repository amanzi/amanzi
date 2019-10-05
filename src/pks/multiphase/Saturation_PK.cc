/*
  This is the multiphase flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
*/

#include <vector>
#include <iostream>

// TPLs
#include "Epetra_Import.h"
#include "Epetra_Vector.h"

// Amanzi
#include "errors.hh"
#include "exceptions.hh"
#include "LinearOperatorFactory.hh"
#include "MFD3D_Diffusion.hh"
#include "OperatorDefs.hh"
#include "TimerManager.hh"

// Amanzi::Multiphase
#include "MultiphaseDefs.hh"
#include "Saturation_PK.hh"

// namespace
namespace Amanzi {
namespace Multiphase {

//Simple construction for multiphase saturation pk
Saturation_PK::Saturation_PK(Teuchos::ParameterList& pk_tree,
                      const Teuchos::RCP<Teuchos::ParameterList>& glist,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& soln):
              bc_saturation(NULL),
              ti_specs_(NULL),
              vo_(NULL),
              src_sink_(NULL),
              S_(S),
              soln_(soln),
              rho_(1.0), mu_(1.0),
              passwd_("state")
{
  linear_operator_list_ = Teuchos::sublist(glist, "Solvers", true);
  preconditioner_list_ = Teuchos::sublist(glist, "Preconditioners", true);
  Teuchos::RCP<Teuchos::ParameterList> flow_list = Teuchos::sublist(glist, "Flow", true);
  mp_list_ = Teuchos::sublist(flow_list, "Saturation Equation", true);
  ti_list_ = Teuchos::sublist(mp_list_, "time integrator", true);
  wrm_list_ = Teuchos::sublist(mp_list_, "Water retention models", true);
  include_capillary_ = flow_list->get<bool>("Include capillary pressure", false);
  jacobian_type_ = flow_list->get<std::string>("Jacobian type", "analytic");
}

void Saturation_PK::Setup()
{
  mesh_ = S_->GetMesh();
  dim_ = mesh_->space_dimension();
  ls_itrs_ = 0;

  // for creating fields
  std::vector<std::string> names(2);
  names[0] = "cell"; 
  names[1] = "face";

  std::vector<AmanziMesh::Entity_kind> locations(2);
  locations[0] = AmanziMesh::CELL; 
  locations[1] = AmanziMesh::FACE;

  std::vector<int> ndofs(2, 1);

  // require state variables for the Multiphase Saturation
  if (!S_->HasField("water_density")) {
    S_->RequireScalar("water_density", passwd_);
  }
  if (!S_->HasField("water_viscosity")) {
    S_->RequireScalar("water_viscosity", passwd_);
  }
  if (!S_->HasField("oil_density")) {
    S_->RequireScalar("oil_density", passwd_);
  }
  if (!S_->HasField("oil_viscosity")) {
    S_->RequireScalar("oil_viscosity", passwd_);
  }
  /*
  if (!S_->HasField("porosity")) {
    S_->RequireScalar("porosity", passwd_);
  }
  */

  if (!S_->HasField("gravity")) {
    S_->RequireConstantVector("gravity", passwd_, dim_);
  }

  if (!S_->HasField("porosity")) {
    S_->RequireField("porosity", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S_->HasField("water_saturation")) {
    S_->RequireField("water_saturation", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S_->HasField("permeability")) {
    S_->RequireField("permeability", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim_);
  }

  if (!S_->HasField("darcy_flux")) {
    S_->RequireField("darcy_flux", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponents(names, locations, ndofs);
//      ->AddComponent("face", AmanziMesh::FACE, 1);
  }
}


/* ******************************************************************
* Clean memory.
****************************************************************** */
Saturation_PK::~Saturation_PK()
{
  if (bc_saturation != NULL) delete bc_saturation;

  if (src_sink_ != NULL) delete src_sink_;
  if (vo_ != NULL) delete vo_;
}


/* ******************************************************************
* 
****************************************************************** */
void Saturation_PK::InitializeFields()
{
  // set popular default values
  if (!S_->GetField("water_density", passwd_)->initialized()) {
    *(S_->GetScalarData("water_density", passwd_)) = 1.0;
    S_->GetField("water_density", passwd_)->set_initialized();
  }

  if (!S_->GetField("water_viscosity", passwd_)->initialized()) {
    *(S_->GetScalarData("water_viscosity", passwd_)) = 1.0;
    S_->GetField("water_viscosity", passwd_)->set_initialized();
  }

  if (!S_->GetField("oil_density", passwd_)->initialized()) {
    *(S_->GetScalarData("oil_density", passwd_)) = 1.0;
    S_->GetField("oil_density", passwd_)->set_initialized();
  }

  if (!S_->GetField("oil_viscosity", passwd_)->initialized()) {
    *(S_->GetScalarData("oil_viscosity", passwd_)) = 1.0;
    S_->GetField("oil_viscosity", passwd_)->set_initialized();
  }

  /*
  if (!S_->GetField("porosity", passwd_)->initialized()) {
    *(S_->GetScalarData("porosity", passwd_)) = 1.0;
    S_->GetField("porosity", passwd_)->set_initialized();
  }
  */

  if (!S_->GetField("gravity", passwd_)->initialized()) {
    Epetra_Vector& gvec = *S_->GetConstantVectorData("gravity", passwd_);
    gvec.PutScalar(0.0);
    if (dim_ == 3) { gvec[dim_ - 1] = -9.80; }
    S_->GetField("gravity", passwd_)->set_initialized();
  }

  if (S_->HasField("porosity")) {
    if (!S_->GetField("porosity", passwd_)->initialized()) {
      S_->GetFieldData("porosity", passwd_)->PutScalar(1.0);
      S_->GetField("porosity", passwd_)->set_initialized();
    }
  }

  if (S_->HasField("water_saturation")) {
    if (!S_->GetField("water_saturation", passwd_)->initialized()) {
      S_->GetFieldData("water_saturation", passwd_)->PutScalar(0.0);
      S_->GetField("water_saturation", passwd_)->set_initialized();
    }
  }

  if (!S_->GetField("permeability", passwd_)->initialized()) {
    S_->GetFieldData("permeability", passwd_)->PutScalar(1.0);
    S_->GetField("permeability", passwd_)->set_initialized();
  }

  if (!S_->GetField("darcy_flux", passwd_)->initialized()) {
    Errors::Message msg;
    msg << "Multiphase Saturation PK: darcy_flux not initialized. You need to initialize it or get it from another PK.\n";
    Exceptions::amanzi_throw(msg);
  } else {
    darcy_flux_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("darcy_flux")));
  }
}


/* ******************************************************************
* 
****************************************************************** */
void Saturation_PK::InitializeSaturation()
{
  // Initilize various common data depending on mesh and state.
  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  K_.resize(ncells_wghost);

  // Fundamental physical quantities
  double* gravity_data;
  S_->GetConstantVectorData("gravity")->ExtractView(&gravity_data);
  gravity_.set(dim_,&(gravity_data[0])); // do it in complicated way because we
                                        // are not sure if gravity_data is an
                                        // array or vector
  g_ = fabs(gravity_[dim_ - 1]);

  // Other constant (temporarily) physical quantaties
  rho_[0] = *(S_->GetScalarData("water_density"));
  rho_[1] = *(S_->GetScalarData("oil_density"));
  mu_[0] = *(S_->GetScalarData("water_viscosity"));
  mu_[1] = *(S_->GetScalarData("oil_viscosity"));
  // phi_ = *S_->GetScalarData("porosity");

  // Initialize miscalleneous default parameters.
  ti_specs_ = NULL;
  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;
  ti_phase_counter = 0;

  src_sink_ = NULL;
  src_sink_distribution_ = 0;

  // create verbosity object. This must precede ComputeBCs(), since
  // vo_ object is used there.
  vo_ = new VerboseObject("MultiphasePK::Saturation", *mp_list_); 

  // Time control specific to this PK.
  ResetPKtimes(0.0, FLOW_INITIAL_DT);
  dT_desirable_ = dT;

  // prepare fractional flow
  rel_perm_w_ = Teuchos::rcp(new RelativePermeability(mesh_));
  rel_perm_n_ = Teuchos::rcp(new RelativePermeability(mesh_));
  rel_perm_w_->Init("wetting", *wrm_list_);
  rel_perm_n_->Init("non wetting", *wrm_list_);

  capillary_pressure_ = Teuchos::rcp(new CapillaryPressure(mesh_));
  capillary_pressure_->Init(*wrm_list_); 

  std::string krel_method_name = mp_list_->get<std::string>("relative permeability");

  // Process parameter list. This entails creating boundary conditions
  // and source term objects.
  ProcessParameterList(*mp_list_);

  // Now compute BCs
  bc_submodel.resize(nfaces_wghost, 0);

  double time = S_->time();
  if (time >= 0.0) T_physics = time;
  time = T_physics;
  bc_saturation->Compute(time, time);
  for (int i = 0; i < bcs_.size(); i++) {
    bcs_[i]->Compute(time, time);
    bcs_[i]->ComputeSubmodel(mesh_);
  }
  ComputeBCs();
  ComputeBC_Pc();

  // create operator boundary condition object
  op_bc_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  op_bc_pc_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  
  // create advection operators
  Teuchos::ParameterList oplist_matrix = mp_list_->sublist("operators").sublist("advection operator");
  op_matrix_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(oplist_matrix, mesh_));
  Teuchos::ParameterList olist_diff = mp_list_->sublist("operators").sublist("diffusion operator").sublist("matrix");
  Operators::PDE_DiffusionFactory opfactory;
  op1_matrix_ = opfactory.Create(olist_diff, mesh_, op_bc_pc_, rho_[1] - rho_[0], gravity_);
  //op1_matrix_->SetDensity(rho_[1] - rho_[0]); // cannot remember why is it minus?
  op_preconditioner_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(oplist_matrix, mesh_));
  op1_preconditioner_ = opfactory.Create(olist_diff, mesh_, op_bc_, 1.0, gravity_);

  // initialize mobility of the wetting phase, total mobility, and fractional flow from relative permeability
  tot_mobility_ = Teuchos::rcp(new CompositeVector(rel_perm_w_->Krel()->Map()));
  tot_mobility_->PutScalar(0.0);

  fractional_flow_ = Teuchos::rcp(new CompositeVector(rel_perm_w_->Krel()->Map()));
  fractional_flow_->PutScalar(0.0);

  dfw_dS_ = Teuchos::rcp(new CompositeVector(rel_perm_w_->Krel()->Map()));
  dfw_dS_->PutScalar(0.0);

  frac_flow_ = 
      Teuchos::rcp(new FractionalFlow(mesh_, mu_[0], rel_perm_w_, mu_[1], rel_perm_n_));
  frac_flow_->Init(*wrm_list_);

  // initialize solution
  sol_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("water_saturation",passwd_)));
  soln_->SetData(sol_);

  // Initalize upwind operator
  Teuchos::ParameterList upw_list = mp_list_->sublist("operators")
                                            .sublist("upwind");
  Operators::UpwindFactory<FractionalFlow> upwind_factory_frac_flow;
  upwind_ = upwind_factory_frac_flow.Create(mesh_, frac_flow_, upw_list);

  Operators::UpwindFactory<RelativePermeability> upwind_factory_rel_perm;
  upwind_n_ = upwind_factory_rel_perm.Create(mesh_, rel_perm_n_, upw_list);

  Operators::UpwindFactory<CapillaryPressure> upwind_pc_factory;
  upwind_pc_ = upwind_pc_factory.Create(mesh_, capillary_pressure_, upw_list);

  // Timer for profiling
  Amanzi::timer_manager.add("ApplyPreconditioner", Amanzi::Timer::ACCUMULATE);
}


void Saturation_PK::InitTimeInterval(Teuchos::ParameterList& ti_list){
  
  ProcessSublistTimeInterval(ti_list, ti_specs_generic_);
 
  ti_specs_generic_.T0  = ti_list.get<double>("start interval time", 0);
  ti_specs_generic_.dT0 = ti_list.get<double>("initial time step", 1);

  double T0 = ti_specs_generic_.T0;
  double dT0 = ti_specs_generic_.dT0;

  dT = dT0;
  dTnext = dT0;

  std::cout<<"T0 "<<T0<<" dT0 "<<dT0<<"\n";

  //if (ti_specs != NULL) OutputTimeHistory(mp_list_, ti_specs->dT_history);
  ti_specs_ = &ti_specs_generic_;

  error_control_ = ti_specs_->error_control_options;
  if (error_control_ == 0) {
    error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE +  // usually 1 [Pa]
                     FLOW_TI_ERROR_CONTROL_SATURATION;  // usually 1e-4;
    ti_specs_->error_control_options = error_control_;
  }

  InitNextTI(T0, dT0, ti_specs_generic_);

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << std::endl 
        << vo_->color("green") << "Initialization Completed!!" << std::endl;
  }
}


/* ******************************************************************
* Generic initialization of a next time integration phase.
****************************************************************** */
void Saturation_PK::InitNextTI(double T0, double dT0, TI_Specs& ti_specs)
{
  std::vector<int>& bc_model = op_bc_->bc_model();
  std::vector<double>& bc_value = op_bc_->bc_value();

  ResetPKtimes(T0, dT0);

  // Compute fractional flow
  Teuchos::RCP<const CompositeVector> water_saturation = S_->GetFieldData("water_saturation");
  frac_flow_->Compute(*water_saturation);
  fractional_flow_ = frac_flow_->Frac_Flow();
  dfw_dS_ = frac_flow_->dF_dS();

  // Now upwind fractional flow
  upwind_->Compute(*darcy_flux_, *darcy_flux_, bc_model, *fractional_flow_);
  upwind_->Compute(*darcy_flux_, *darcy_flux_, bc_model, *dfw_dS_); 

  fractional_flow_->Multiply(1.0, *fractional_flow_, *darcy_flux_, 0.0);
  dfw_dS_->Multiply(1.0, *dfw_dS_, *darcy_flux_, 0.0);

  // Initialize operators
  op_matrix_->Setup(*darcy_flux_);
  op_matrix_->UpdateMatrices(fractional_flow_.ptr(), Teuchos::null);

  op_preconditioner_->Setup(*darcy_flux_);
  op_preconditioner_->UpdateMatrices(dfw_dS_.ptr(), Teuchos::null);
  // op_preconditioner_->UpdateMatrices(*fractional_flow_);
  // op_preconditioner_->ApplyBCs(op_bc_, true);
  // op_preconditioner_->global_operator()->SymbolicAssembleMatrix();
  // op_ = op_preconditioner_->global_operator();
  // op_->SymbolicAssembleMatrix();

  SetAbsolutePermeabilityTensor();
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
  op1_matrix_->Setup(Kptr, Teuchos::null, Teuchos::null);
  op1_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op1_matrix_->ApplyBCs(true, true, true);

  op1_preconditioner_->Setup(Kptr, fractional_flow_, Teuchos::null);
  op1_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // Initialize source
  // if (src_sink_ != NULL) {
  //   src_sink_->ComputeDistribute(T0, T0 + dT0);
  // }

  // set up new time integration or solver
  std::string ti_method_name(ti_specs.ti_method_name);

  if (ti_specs.ti_method == FLOW_TIME_INTEGRATION_BDF1) {
    Teuchos::ParameterList bdf1_list = ti_specs.ti_list_ptr_->sublist("BDF1");
    if (! bdf1_list.isSublist("VerboseObject"))
        bdf1_list.sublist("VerboseObject") = mp_list_->sublist("VerboseObject");

    bdf1_dae = Teuchos::rcp(new BDF1_TI<TreeVector, TreeVectorSpace>(*this, bdf1_list, soln_));
  }
}


/* ******************************************************************* 
* Performs one time integration step of size dT. The boundary conditions 
* and source terms are calculated only once at the beginning.  
******************************************************************* */
bool Saturation_PK::Advance(double dT_MPC, double& dT_actual)
{
  dT = dT_MPC;
  double time = S_->time();
  if (time >= 0.0) T_physics = time;
  time = T_physics;

  // Update boundary conditions and source terms for each time iteration
  // if (src_sink_ != NULL) {
  //   src_sink_->ComputeDistribute(time, time + dT_MPC);
  // }

  bc_saturation->Compute(time, time);
  for (int i = 0; i < bcs_.size(); i++) {
    bcs_[i]->Compute(time, time);
    bcs_[i]->ComputeSubmodel(mesh_);
  }
  ComputeBCs();

  // Update the flux
  darcy_flux_ = S_->GetFieldData("darcy_flux", passwd_);

  Teuchos::RCP<TreeVector> udot = Teuchos::rcp(new TreeVector(soln_->Map()));
  udot->PutScalar(0.0);

  // predict water mass change during time step
  // if (ti_specs_->num_itrs == 0) {  // initialization
    // I do not know how to estimate du/dt, so it is set to zero.
    // Teuchos::RCP<TreeVector> udot = Teuchos::rcp(new TreeVector(soln_->Map()));
    // udot->PutScalar(0.0);
    /*
    if (ti_specs_->ti_method == FLOW_TIME_INTEGRATION_BDF1) {
      std::cout << "Set Initial State in Advance, soln_ = " << *soln_->Data()->ViewComponent("cell") << "\n";
      bdf1_dae->SetInitialState(time, soln_, udot);
    } else {
      Errors::Message msg;
      msg << "Multiphase Saturation PK: Time integration only supports BDF1 at the moment \n";
      Exceptions::amanzi_throw(msg);
    }
    */
    //UpdatePreconditioner(time, soln_, dT);
    //ti_specs_->num_itrs++;
  // }

  /*
  if (ti_specs_->ti_method == FLOW_TIME_INTEGRATION_BDF1) {
    while (bdf1_dae->TimeStep(dT, dTnext, soln_)) {
      dT = dTnext;
    }
    // --etc this is a bug in general -- should not commit the solution unless
    //   the step passes all other PKs
    bdf1_dae->CommitSolution(dT, soln_);
    T_physics = bdf1_dae->time();
  }
  */

  bdf1_dae->SetInitialState(time, soln_, udot);

  bool fail = false;
  if (ti_specs_->ti_method == FLOW_TIME_INTEGRATION_BDF1){
    fail = bdf1_dae->TimeStep(dT, dTnext, soln_);
    if (fail) {
      dT = dTnext;
      return fail;
    }
  }

  // tell the caller what time step we actually took
  dT_actual = dT;
  
  dt_tuple times(time, dT);
  ti_specs_->dT_history.push_back(times);

  ti_specs_->num_itrs++;

  dT = dTnext;

  return fail;
}


/* ******************************************************************* 
* Wrapper for Saturation_PK::Advance. This is the interface used by
* the MPC.
******************************************************************* */
bool Saturation_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool failed = false;
  double dt = t_new - t_old;
  double dt_actual(dt);
  failed = this->Advance(dt, dt_actual);
  return failed;
}


/* ******************************************************************* 
* Write back the solution and other data to the State 
******************************************************************* */
void Saturation_PK::CommitState(const Teuchos::Ptr<State>& S)
{
  Epetra_MultiVector& ws = *S->GetFieldData("water_saturation", passwd_)->ViewComponent("cell");
  Epetra_MultiVector& sol_c = *soln_->Data()->ViewComponent("cell"); 
  ws = sol_c;
  //std::cout << "solution: " << ws << "\n";
}


/* ******************************************************************* 
* Wrapper for Saturation_PK::CommitState. This is the interface used 
* by the MPC.
******************************************************************* */
void Saturation_PK::CommitStep(double t_old, double t_new)
{
  this->CommitState(S_.ptr());
}


/* ******************************************************************* 
* Process the source term particular for the saturation equation.
* This is a simple model for wells. Injection wells are assumed to 
* consist of only one phase (e.g. water). Production wells only produce 
* the other phase (e.g. oil).
******************************************************************* */
void Saturation_PK::AddSourceTerm(CompositeVector& rhs, const CompositeVector& fractional_flow)
{
  /*
  Epetra_MultiVector& rhs_cell = *rhs.ViewComponent("cell");
  const Epetra_MultiVector& fractional_flow_c = *fractional_flow.ViewComponent("cell");
  Flow::FlowDomainFunction::Iterator src;
  for (src = src_sink_->begin(); src != src_sink_->end(); ++src) {
    int c = src->first;
    rhs_cell[0][c] += mesh_->cell_volume(c) * (std::max(src->second, 0.0)
                   + std::min(src->second, 0.0) * fractional_flow_c[0][c]);
  }
  */
}


/* *********************************************************************
* Compute boundary conditions
********************************************************************* */
void Saturation_PK::ComputeBCs()
{
  std::vector<int>& bc_model = op_bc_->bc_model();
  std::vector<double>& bc_value = op_bc_->bc_value();
  std::vector<double>& bc_mixed = op_bc_->bc_mixed();

  std::vector<int>& bc_model_pc = op_bc_pc_->bc_model();
  std::vector<double>& bc_value_pc = op_bc_pc_->bc_value();

  int flag_essential_bc = 0;
  dirichlet_bc_faces_ = 0;

  for (int n = 0; n < bc_model.size(); n++) {
    bc_model[n] = Operators::OPERATOR_BC_NONE;
    bc_value[n] = 0.0;
  }

  Flow::FlowBoundaryFunction::Iterator bc;
  for (bc = bc_saturation->begin(); bc != bc_saturation->end(); ++bc) {
    int f = bc->first;
    bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
    bc_value[f] = bc->second[0];
    bc_model_pc[f] = Operators::OPERATOR_BC_DIRICHLET;
    flag_essential_bc = 1;
    dirichlet_bc_faces_++;
  }

  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->bc_name() == "flux") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value[f] = it->second[0];
        bc_model_pc[f] = Operators::OPERATOR_BC_NEUMANN;
      }
    }
  }

  // mark missing boundary conditions as zero flux conditions 
  AmanziMesh::Entity_ID_List cells;
  missed_bc_faces_ = 0;
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_model[f] == Operators::OPERATOR_BC_NONE) {
      cells.clear();
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int ncells = cells.size();

      if (ncells == 1) {
        bc_model_pc[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value_pc[f] = 0.0;
        missed_bc_faces_++;
      }
    }
  }

  // verify that the algebraic problem is consistent
  #ifdef HAVE_MPI
  int flag = flag_essential_bc;
  mesh_->get_comm()->MaxAll(&flag, &flag_essential_bc, 1);  // find the global maximum
  #endif
  if (! flag_essential_bc && vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "WARNING: no essential boundary conditions, solver may fail" << std::endl;
  }
}


/* ******************************************************************
* Compute boundary conditions for capillary pressure                                          
****************************************************************** */
void Saturation_PK::ComputeBC_Pc()
{
  std::vector<int>& bc_model = op_bc_pc_->bc_model();
  std::vector<double>& bc_value = op_bc_pc_->bc_value();
  std::vector<double>& bc_mixed = op_bc_pc_->bc_mixed();

  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm_n_->WRM(); 

  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED, &block);
    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      if (bc_model[*i] == Operators::OPERATOR_BC_DIRICHLET) {
        bc_value[*i] = WRM[mb]->capillaryPressure(bc_value[*i]);
      } else if (bc_model[*i] == Operators::OPERATOR_BC_NEUMANN) {
        bc_value[*i] = bc_value[*i];
      }
    }
  }
}


/* ******************************************************************
*
****************************************************************** */
void Saturation_PK::SetAbsolutePermeabilityTensor()
{
  const CompositeVector& cv = *S_->GetFieldData("permeability");
  cv.ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& perm = *cv.ViewComponent("cell", true);

  if (dim_ == 2) {
    for (int c = 0; c < K_.size(); c++) {
      if (perm[0][c] == perm[1][c]) {
        K_[c].Init(dim_, 1);
        K_[c](0, 0) = perm[0][c];
      } else {
        K_[c].Init(dim_, 2);
        K_[c](0, 0) = perm[0][c];
        K_[c](1, 1) = perm[1][c];
      }
    }    
  } else if (dim_ == 3) {
    for (int c = 0; c < K_.size(); c++) {
      if (perm[0][c] == perm[1][c] && perm[0][c] == perm[2][c]) {
        K_[c].Init(dim_, 1);
        K_[c](0, 0) = perm[0][c];
      } else {
        K_[c].Init(dim_, 2);
        K_[c](0, 0) = perm[0][c];
        K_[c](1, 1) = perm[1][c];
        K_[c](2, 2) = perm[2][c];
      }
    }        
  }
}

}  // namespace Multiphase
}  // namespace Amanzi
