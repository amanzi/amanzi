/*
  Multiphase

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
#include "CommonDefs.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "LinearOperatorFactory.hh"
#include "Mesh.hh"
#include "MFD3D_Diffusion.hh"
#include "OperatorDefs.hh"
#include "PK_DomainFunctionFactory.hh"

// Amanzi::Multiphase
#include "CompH_PK.hh"
#include "Multiphase_Utils.hh"
#include "MultiphaseDefs.hh"

// namespace
namespace Amanzi {
namespace Multiphase {

CompH_PK::CompH_PK(Teuchos::ParameterList& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& specs_list,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& soln)
  : src_sink_(NULL),
    S_(S),
    soln_(soln),
    passwd_("state"),
    jacobian_type_("analytic")
{
  comp_list_ = Teuchos::sublist(specs_list, "Component 2");
  if (!specs_list->isSublist("Water retention models")) {
    Errors::Message msg;
    msg << "Component PK: input XML files does not have <Water retention models> sublist in <MPMC Specs> \n";
    Exceptions::amanzi_throw(msg);
  } else {
    wrm_list_ = specs_list->sublist("Water retention models");  
  }

  jacobian_type_ = specs_list->get<std::string>("Jacobian type", "analytic");
  
  Cg_ = specs_list->sublist("Component 2").get<double>("compressibility constant");
  /*
  if (!comp_list_->isParameter("component name")) {
    Errors::Message msg("CompH_PK: parameter \"component name\" is missing");
    Exceptions::amanzi_throw(msg);
  } else {
    comp_name_ = comp_list_->get<int>("component id");   
  }
  if (!comp_list_->isParameter("component coefficient")) {
    Errors::Message msg("CompH_PK: parameter \"component coefficient\" is missing");
    Exceptions::amanzi_throw(msg);
  } else {
    comp_coeff_ = comp_list_->get<double>("component coefficient");
  }
  */

  mesh_ = S->GetMesh();
  dim_ = mesh_->space_dimension();
  src_sink_ = NULL;
  src_sink_distribution_ = 0;
    
  /*
  std::vector<std::string> names(2);
  names[0] = "cell";
  names[1] = "face";

  std::vector<int> ndofs(2, 1);

  std::vector<AmanziMesh::Entity_kind> locations(2);
  locations[0] = AmanziMesh::CELL;
  locations[1] = AmanziMesh::FACE;
  */

  // require state variables for the MultiphasePhase1 
  if (!S->HasField("phase1_density")) {
    S->RequireScalar("phase1_density", passwd_);
  }

  if (!S->HasField("phase1_viscosity")) {
    S->RequireScalar("phase1_viscosity", passwd_);
  }

  if (!S->HasField("phase2_viscosity")) {
    S->RequireScalar("phase2_viscosity", passwd_);
  }

  if (!S->HasField("porosity")) {
    S->RequireScalar("porosity", passwd_);
  }

  if (!S->HasField("gravity")) {
    S->RequireConstantVector("gravity", passwd_, dim_);
  }

  if (!S->HasField("pressure_w")) {
    S->RequireField("pressure_w", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  
  if (!S->HasField("pressure_n")) {
    S->RequireField("pressure_n", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField("permeability")) {
    S->RequireField("permeability", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim_);
  }

  if (!S->HasField("diffusion tensor")) {
    S->RequireField("diffusion tensor", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim_);
  }

  if (!S->HasField("saturation_w")) {
    S->RequireField("saturation_w", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField("saturation_n")) {
    S->RequireField("saturation_n", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField("hydrogen density liquid")) {
    S->RequireField("hydrogen density liquid", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField("dPc_dS")) {
    S->RequireField("dPc_dS", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField("velocity_wet")) {
    S->RequireField("velocity_wet", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      //->SetComponents(names, locations, ndofs);
      ->SetComponent("face", AmanziMesh::FACE, 1);
  }

  if (!S->HasField("velocity_nonwet")) {
    S->RequireField("velocity_nonwet", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);
  }
}


/* ******************************************************************
* Clean memory.
****************************************************************** */
CompH_PK::~CompH_PK()
{
  if (src_sink_ != NULL) delete src_sink_;
  if (vo_ != NULL) delete vo_;
}


/* ******************************************************************
* Set popular default values
****************************************************************** */
void CompH_PK::InitializeFields()
{
  if (!S_->GetField("phase1_viscosity", passwd_)->initialized()) {
    *(S_->GetScalarData("phase1_viscosity", passwd_)) = 1.0;
    S_->GetField("phase1_viscosity", passwd_)->set_initialized();
  }

  if (!S_->GetField("phase2_viscosity", passwd_)->initialized()) {
    *(S_->GetScalarData("phase2_viscosity", passwd_)) = 1.0;
    S_->GetField("phase2_viscosity", passwd_)->set_initialized();
  }

  if (!S_->GetField("porosity", passwd_)->initialized()) {
    *(S_->GetScalarData("porosity", passwd_)) = 1.0;
    S_->GetField("porosity", passwd_)->set_initialized();
  }

  if (!S_->GetField("gravity", passwd_)->initialized()) {
    Epetra_Vector& gvec = *S_->GetConstantVectorData("gravity", passwd_);
    gvec.PutScalar(0.0);
    if (dim_ == 3) { gvec[dim_ - 1] = -9.80; }
    S_->GetField("gravity", passwd_)->set_initialized();
  }

  if (!S_->GetField("permeability", passwd_)->initialized()) {
    S_->GetFieldData("permeability", passwd_)->PutScalar(1.0);
    S_->GetField("permeability", passwd_)->set_initialized();
  }

  if (!S_->GetField("diffusion tensor", passwd_)->initialized()) {
    S_->GetFieldData("diffusion tensor", passwd_)->PutScalar(1.0);
    S_->GetField("diffusion tensor", passwd_)->set_initialized();
  }

  if (S_->HasField("saturation_n")) {
    if (!S_->GetField("saturation_n", passwd_)->initialized()) {
      S_->GetFieldData("saturation_n", passwd_)->PutScalar(0.0);
      S_->GetField("saturation_n", passwd_)->set_initialized();
    }
  }

  if (!S_->GetField("phase1_density", passwd_)->initialized()) {
    *(S_->GetScalarData("phase1_density", passwd_)) = 1.0;
    S_->GetField("phase1_density", passwd_)->set_initialized();
  }

  if (S_->HasField("hydrogen density liquid")) {
    if (!S_->GetField("hydrogen density liquid", passwd_)->initialized()) {
      S_->GetFieldData("hydrogen density liquid", passwd_)->PutScalar(0.0);
      S_->GetField("hydrogen density liquid", passwd_)->set_initialized();
    }
  }

  if (!S_->GetField("pressure_w", passwd_)->initialized()) {
    S_->GetFieldData("pressure_w", passwd_)->PutScalar(0.0);
    S_->GetField("pressure_w", passwd_)->set_initialized();
  }

  if (!S_->GetField("pressure_n", passwd_)->initialized()) {
    S_->GetFieldData("pressure_n", passwd_)->PutScalar(0.0);
    S_->GetField("pressure_n", passwd_)->set_initialized();
  }

  if (!S_->GetField("dPc_dS", passwd_)->initialized()) {
    S_->GetFieldData("dPc_dS", passwd_)->PutScalar(0.0);
    S_->GetField("dPc_dS", passwd_)->set_initialized();
  }

  if (!S_->GetField("velocity_wet", passwd_)->initialized()) {
    S_->GetFieldData("velocity_wet", passwd_)->PutScalar(0.0);
    S_->GetField("velocity_wet", passwd_)->set_initialized();
  }

  if (!S_->GetField("velocity_nonwet", passwd_)->initialized()) {
    S_->GetFieldData("velocity_nonwet", passwd_)->PutScalar(0.0);
    S_->GetField("velocity_nonwet", passwd_)->set_initialized();
  } 
}


void CompH_PK::InitializeComponent()
{
  // Initilize various common data depending on mesh and state.
  ncells_owned_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncells_wghost_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  nfaces_owned_ = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nfaces_wghost_ = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

  // Fundamental physical quantities
  double* gravity_data;
  S_->GetConstantVectorData("gravity")->ExtractView(&gravity_data);
  gravity_.set(dim_,&(gravity_data[0])); // do it in complicated way because we
                                        // are not sure if gravity_data is an
                                        // array or vector
  g_ = fabs(gravity_[dim_ - 1]);

  // Other constant (temporarily) physical quantaties
  mu1_ = *S_->GetScalarData("phase1_viscosity");
  mu2_ = *S_->GetScalarData("phase2_viscosity");
  phi_ = *S_->GetScalarData("porosity");
  rho_ = *S_->GetScalarData("phase1_density");
  rho_w_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("water component density")));

  // Initialize miscalleneous default parameters.
  ti_specs_ = NULL;
  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;

  src_sink_ = NULL;
  src_sink_distribution_ = 0;

  T_physics_ = 0.0;

  // create verbosity object
  vo_ = new VerboseObject("MPMC:Component", *comp_list_); 

  // Allocate memory for boundary data.
  Teuchos::RCP<Flow::FlowBoundaryFunction> bc;
  auto& bc_list = comp_list_->sublist("boundary conditions");

  bcs_.clear();

  // -- pressure 
  if (bc_list.isSublist("pressure")) {
    PK_DomainFunctionFactory<Flow::FlowBoundaryFunction> bc_factory(mesh_);

    Teuchos::ParameterList& tmp_list = bc_list.sublist("pressure");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);
        bc = bc_factory.Create(spec, "boundary pressure", AmanziMesh::FACE, Teuchos::null);
        bc->set_bc_name("pressure");
        bcs_.push_back(bc);
      }
    }
  }

  // -- saturation
  if (bc_list.isSublist("saturation")) {
    PK_DomainFunctionFactory<Flow::FlowBoundaryFunction> bc_factory(mesh_);

    Teuchos::ParameterList& tmp_list = bc_list.sublist("saturation");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);
        bc = bc_factory.Create(spec, "boundary saturation", AmanziMesh::FACE, Teuchos::null);
        bc->set_bc_name("saturation");
        bcs_.push_back(bc);
      }
    }
  }

  // -- hydrogen density
  if (bc_list.isSublist("hydrogen density")) {
    PK_DomainFunctionFactory<Flow::FlowBoundaryFunction> bc_factory(mesh_);

    Teuchos::ParameterList& tmp_list = bc_list.sublist("hydrogen density");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);
        bc = bc_factory.Create(spec, "boundary hydrogen density", AmanziMesh::FACE, Teuchos::null);
        bc->set_bc_name("hydrogen density");
        bcs_.push_back(bc);
      }
    }
  }

  // -- mass flux
  if (bc_list.isSublist("mass flux total")) {
    PK_DomainFunctionFactory<Flow::FlowBoundaryFunction> bc_factory(mesh_);

    Teuchos::ParameterList& tmp_list = bc_list.sublist("mass flux total");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      if (it->second.isList()) {
        Teuchos::ParameterList spec = Teuchos::getValue<Teuchos::ParameterList>(it->second);
        bc = bc_factory.Create(spec, "outward mass flux", AmanziMesh::FACE, Teuchos::null);
        bc->set_bc_name("flux");
        bcs_.push_back(bc);
      }
    }
  }

  // double time = S_->time();
  double time = T_physics_;
  for (int i = 0; i < bcs_.size(); i++) {
    bcs_[i]->Compute(time, time);
  }

  // Create water retention models.
  coef_w_ = Teuchos::rcp(new MPCoeff(mesh_));
  coef_n_ = Teuchos::rcp(new MPCoeff(mesh_));
  coef_w_->Init("wetting", wrm_list_);
  coef_n_->Init("non wetting", wrm_list_, Cg_);
  capillary_pressure_ = Teuchos::rcp(new CapillaryPressure(mesh_));
  capillary_pressure_->Init(wrm_list_); 

  // create and ppulate boundary conditions
  op_bc_p_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  op_bc_p_n_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  op_bc_rhl_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  op_bc_s_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

  // allocate memory for absolute permeability
  K_.resize(ncells_wghost_);
  D1_.resize(ncells_wghost_);
  ConvertFieldToTensor(S_, dim_, "permeability", K_);
  ConvertFieldToTensor(S_, dim_, "diffusion tensor", D1_);

  // Select a proper matrix class. 
  Teuchos::ParameterList& tmp_list = comp_list_->sublist("operators").sublist("diffusion operator");
  Teuchos::ParameterList oplist_matrix = tmp_list.sublist("matrix");
  Teuchos::ParameterList oplist_pc = tmp_list.sublist("preconditioner");
  op_list_ = Teuchos::rcp(new Teuchos::ParameterList(tmp_list.sublist("preconditioner")));

  // create diffusion operators op_matrix_ and op_preconditioner_.
  // They will need to be initialized, which is done later in InitNextTI.
  Operators::PDE_DiffusionFactory opfactory;
  op1_matrix_ = Teuchos::rcp(new Operators::PDE_DiffusionFVwithGravity(oplist_matrix, mesh_, 1.0, gravity_));
  op1_matrix_->SetBCs(op_bc_p_, op_bc_p_);

  op2_matrix_ = Teuchos::rcp(new Operators::PDE_DiffusionFVwithGravity(oplist_matrix, mesh_, 1.0, gravity_));
  op2_matrix_->SetBCs(op_bc_p_n_, op_bc_p_n_);

  op3_matrix_ = opfactory.Create(oplist_matrix, mesh_, op_bc_rhl_, 1.0, gravity_);

  ComputeBCs(false);
  ComputeBC_Pn();

  // preconditioner operators
  op1_preconditioner_ = opfactory.Create(oplist_pc, mesh_, op_bc_p_, 1.0, gravity_);
  op2_preconditioner_ = opfactory.Create(oplist_pc, mesh_, op_bc_rhl_, 1.0, gravity_);
  op3_preconditioner_ = opfactory.Create(oplist_pc, mesh_, op_bc_s_, 1.0, gravity_);

  Teuchos::ParameterList olist_adv = comp_list_->sublist("operators").sublist("advection operator");
  op_prec_pres_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, op1_preconditioner_->global_operator()));
  op_prec_sat_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, op3_preconditioner_->global_operator()));
  op_prec_rho_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, op2_preconditioner_->global_operator()));
  op_prec_pres_->SetBCs(op_bc_p_n_, op_bc_p_n_);
  op_prec_sat_->SetBCs(op_bc_s_, op_bc_s_);
  op_prec_rho_->SetBCs(op_bc_p_, op_bc_p_);

  // accumulation terms
  op_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op_prec_pres_->global_operator()));
  op1_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op_prec_rho_->global_operator()));
  op2_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op_prec_sat_->global_operator()));

  if (jacobian_type_ == "numerical")
    for (int n_op = 0; n_op < soln_->size(); n_op++) {
      ops_.push_back(opfactory.Create(oplist_pc, mesh_, op_bc_p_, 1.0, gravity_));
    }

  // init upwind operator
  Teuchos::ParameterList upw_list = comp_list_->sublist("operators")
                                               .sublist("diffusion operator")
                                               .sublist("upwind");

  upwind_n_ = Teuchos::rcp(new Operators::UpwindFlux<MPCoeff>(mesh_, coef_n_));
  upwind_n_->Init(upw_list);
  upwind_w_ = Teuchos::rcp(new Operators::UpwindFlux<MPCoeff>(mesh_, coef_w_));
  upwind_w_->Init(upw_list);

  // init upwind_velocity
  upwind_vw_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("velocity_wet", passwd_)));
  upwind_vn_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("velocity_nonwet", passwd_)));

  tmp_flux1_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("velocity_wet", passwd_)));
  tmp_flux2_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("velocity_nonwet", passwd_)));
}


/* ******************************************************************
* Generic initialization of a next time integration phase.
****************************************************************** */
void CompH_PK::InitNextTI()
{
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > D1ptr = Teuchos::rcpFromRef(D1_);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->Setup(Kptr, coef_w_->Krel(), Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);

  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->Setup(Kptr, coef_n_->Krel(), Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);

  op3_matrix_->Setup(D1ptr, Teuchos::null, Teuchos::null);
  op3_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);

  if (jacobian_type_ == "numerical")
    for (int ii = 0; ii < ops_.size(); ii++) {
      ops_[ii]->Setup(Kptr, coef_n_->Krel(), Teuchos::null);
      ops_[ii]->UpdateMatrices(Teuchos::null, Teuchos::null);
    }
} // End InitNextTI


/* ******************************************************************* 
* Wrapper for CompH_PK::CommitState. This is the interface used 
* by the MPC.
******************************************************************* */
void CompH_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  const std::vector<int>& bc_model_p_n = op_bc_p_n_->bc_model();
  const std::vector<double>& bc_value_p_n = op_bc_p_n_->bc_value();

  // Write pressure to state
  auto p2 = S_->GetFieldData("pressure_n", passwd_);
  CompositeVector& p1 = *S_->GetFieldData("pressure_w", passwd_);
  p1 = *soln_->SubVector(0)->Data();

  // write saturation of phase1 to state
  CompositeVector& s1 = *S_->GetFieldData("saturation_w", passwd_);
  s1 = *soln_->SubVector(1)->Data();
  CompositeVector& s2 = *S_->GetFieldData("saturation_n", passwd_);
  s2.Update(-1.0, s1, 0.0);
  s2.Shift(1.0);

  capillary_pressure_->Compute(s1); 
  p2->Update(1.0, p1, 0.0);
  p2->Update(1.0, *capillary_pressure_->Pc(), 1.0);

  capillary_pressure_->dPc_dS()->Scale(-1.0);
  CompositeVector& dPc_dS = *S_->GetFieldData("dPc_dS", passwd_);
  dPc_dS = *capillary_pressure_->dPc_dS();

  // compute upwind velocities from pressure 
  // Calculate total mobility needed to initialize diffusion operator
  coef_n_->Compute(s1, bc_model_p_n, bc_value_p_n);
  upwind_vn_ = S_->GetFieldData("velocity_nonwet", passwd_);

  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_p_n, *coef_n_->Krel());
  coef_n_->Krel()->Scale((t_new - t_old)/mu2_); // don't need rho to compute darcy flux

  auto phase2_flux = S_->GetFieldData("velocity_nonwet", passwd_);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->Setup(Kptr, coef_n_->Krel(), Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->ApplyBCs(true, true, true);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->UpdateFlux(p2.ptr(), phase2_flux.ptr());

  // compute upwind velocities from pressure 
  // CompositeVector& phase2_flux = *S_->GetFieldData("velocity_nonwet", passwd_);
  // ((Teuchos::RCP<Operators::PDE_Diffusion>)op2_matrix_)->UpdateFlux(p2, phase2_flux);

  // write hydrogen density in liquid to state
  CompositeVector& rhl = *S_->GetFieldData("hydrogen density liquid", passwd_);
  rhl = *soln_->SubVector(2)->Data();
}


// need to fix flag essential bc. We need Dirichlet bc for all 4 primary
// variables before incrementing Dirichlet bc counter 
/* ******************************************************************* 
* Compute boundary conditions based on bc type and values
******************************************************************* */
void CompH_PK::ComputeBCs(bool stop)
{
  std::vector<int>& bc_model_p = op_bc_p_->bc_model();
  std::vector<double>& bc_value_p = op_bc_p_->bc_value();

  std::vector<int>& bc_model_p_n = op_bc_p_n_->bc_model();
  std::vector<double>& bc_value_p_n = op_bc_p_n_->bc_value();

  std::vector<int>& bc_model_s = op_bc_s_->bc_model();
  std::vector<double>& bc_value_s = op_bc_s_->bc_value();

  std::vector<int>& bc_model_rhl = op_bc_rhl_->bc_model();
  std::vector<double>& bc_value_rhl = op_bc_rhl_->bc_value();

  int flag_essential_bc = 0;
  dirichlet_bc_faces_ = 0;

  for (int n = 0; n < bc_model_p.size(); n++) {
    bc_model_p[n] = Operators::OPERATOR_BC_NONE;

    bc_value_p[n] = 0.0;
    bc_value_p_n[n] = 0.0;
    bc_value_s[n] = 0.0;
    bc_value_rhl[n] = 0.0;
  }

  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->bc_name() == "pressure") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        bc_model_p[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value_p[f] = it->second[0];
        flag_essential_bc = 1;
        dirichlet_bc_faces_++;
      }
    }

    if (bcs_[i]->bc_name() == "saturation") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        bc_value_s[f] = it->second[0];
      }
    }

    if (bcs_[i]->bc_name() == "hydrogen density") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        bc_value_rhl[f] = it->second[0];
      }
    }

    if (bcs_[i]->bc_name() == "flux") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        bc_model_p[f] = Operators::OPERATOR_BC_NEUMANN;
        if (stop) {
          bc_value_p[f] = 0.0;
        } else {
          bc_value_p[f] = it->second[0];
        }
      }
    }
  }

  // mark missing boundary conditions as zero flux conditions 
  AmanziMesh::Entity_ID_List cells;
  missed_bc_faces_ = 0;
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();

    if (ncells == 1) {
      if (bc_model_p[f] == Operators::OPERATOR_BC_NONE) {
        bc_model_p[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value_p[f] = 0.0;
        bc_value_p_n[f] = 0.0;
        bc_value_s[f] = 0.0;
        bc_value_rhl[f] = 0.0;
      }

      missed_bc_faces_++;
    }
  }

  bc_model_p_n = bc_model_p;
  bc_model_s = bc_model_p;
  bc_model_rhl = bc_model_p;

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
void CompH_PK::ComputeBC_Pn()
{
  const std::vector<int>& bc_model_p = op_bc_p_->bc_model();
  const std::vector<double>& bc_value_p = op_bc_p_->bc_value();
  const std::vector<double>& bc_value_s = op_bc_s_->bc_value();

  std::vector<double>& bc_value_p_n = op_bc_p_n_->bc_value();

  AmanziMesh::Entity_ID_List cells;
  auto& wrm = coef_n_->wrm();

  for (int f = 0; f != nfaces_owned_; ++f) {
    if (bc_model_p[f] == Operators::OPERATOR_BC_DIRICHLET) {
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int c = cells[0];

      int mb = (*wrm->first)[c];
      bc_value_p_n[f] = bc_value_p[f] + (wrm->second)[mb]->capillaryPressure(bc_value_s[f]);
    }
  }
}

}  // namespace Flow
} // End namespace Amanzi
