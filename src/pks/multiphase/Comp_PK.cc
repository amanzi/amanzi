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

#include "Epetra_Import.h"
#include "Epetra_Vector.h"

#include "errors.hh"
#include "exceptions.hh"
#include "LinearOperatorFactory.hh"
#include "mfd3d_diffusion.hh"
#include "OperatorDefs.hh"
#include "CommonDefs.hh"

#include "Mesh.hh"

#include "Comp_PK.hh"
#include "MultiphaseDefs.hh"

// namespace
namespace Amanzi {
namespace Multiphase {

Comp_PK::Comp_PK(Teuchos::ParameterList& pk_tree,
                      const Teuchos::RCP<Teuchos::ParameterList>& specs_list,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& soln):
              bc_pressure_(NULL),
              bc_mfw_(NULL),
              bc_mfn_(NULL),
              bc_saturation_(NULL),
              src_sink_(NULL),
              S_(S),
              soln_(soln),
              passwd_("state")
{
  comp_name_ = pk_tree.get<std::string>("component name");
  //std::cout << "component label: " << comp_name_ << "\n";
  comp_list_ = Teuchos::rcp(new Teuchos::ParameterList(specs_list->sublist(comp_name_)));
  molar_fraction_name1_ = "molar_fraction_" + comp_list_->get<std::string>("component name") + "_phase1";
  molar_fraction_name2_ = "molar_fraction_" + comp_list_->get<std::string>("component name") + "_phase2";
  diff_tensor_name1_ = "diff_tensor_" + comp_list_->get<std::string>("component name") + "_phase1";
  diff_tensor_name2_ = "diff_tensor_" + comp_list_->get<std::string>("component name") + "_phase2";

  if (!specs_list->isSublist("Water retention models")) {
    Errors::Message msg;
    msg << "Component PK: input XML files does not have <Water retention models> sublist in <MPMC Specs> \n";
    Exceptions::amanzi_throw(msg);
  } else {
    wrm_list_ = specs_list->sublist("Water retention models");  
  }
  
  if (!comp_list_->isParameter("component id")) {
    Errors::Message msg("Comp_PK: parameter \"component id\" is missing");
    Exceptions::amanzi_throw(msg);
  } else {
    comp_id_ = comp_list_->get<int>("component id");   
  }
  comp_coeff_ = specs_list->sublist("EOS").sublist(comp_name_).get<double>("coefficient value");
  /*
  if (!comp_list_->isParameter("component name")) {
    Errors::Message msg("Comp_PK: parameter \"component name\" is missing");
    Exceptions::amanzi_throw(msg);
  } else {
    comp_name_ = comp_list_->get<int>("component id");   
  }
  if (!comp_list_->isParameter("component coefficient")) {
    Errors::Message msg("Comp_PK: parameter \"component coefficient\" is missing");
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

  if (!S->HasField("phase2_density")) {
    S->RequireScalar("phase2_density", passwd_);
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

  if (!S->HasField(diff_tensor_name1_)) {
    S->RequireField(diff_tensor_name1_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim_);
  }

  if (!S->HasField(diff_tensor_name2_)) {
    S->RequireField(diff_tensor_name2_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
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

  if (!S->HasField(molar_fraction_name1_)) {
    S->RequireField(molar_fraction_name1_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField(molar_fraction_name2_)) {
    S->RequireField(molar_fraction_name2_, passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

 if (!S->HasField("fugacity1")) {
    S->RequireField("fugacity1", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

 if (!S->HasField("fugacity2")) {
    S->RequireField("fugacity2", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField("velocity_wet")) {
    S->RequireField("velocity_wet", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      //->SetComponents(names, locations, ndofs);
      ->SetComponent("face", AmanziMesh::FACE, 1);
  }

  if (!S->HasField("velocity_nonwet")) {
    S->RequireField("velocity_nonwet", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      //->SetComponents(names, locations, ndofs);
      ->SetComponent("face", AmanziMesh::FACE, 1);
  }
} // End Comp_PK constructor


/* ******************************************************************
* Clean memory.
****************************************************************** */
Comp_PK::~Comp_PK()
{
  if (bc_pressure_ != NULL) delete bc_pressure_;
  if (bc_saturation_ != NULL) delete bc_saturation_;
  if (bc_mfw_ != NULL) delete bc_mfw_;
  if (bc_mfn_ != NULL) delete bc_mfn_;

  if (src_sink_ != NULL) delete src_sink_;
  if (vo_ != NULL) delete vo_;
}


void Comp_PK::InitializeFields()
{
  // set popular default values

  if (!S_->GetField("phase1_density", passwd_)->initialized()) {
    *(S_->GetScalarData("phase1_density", passwd_)) = 1.0;
    S_->GetField("phase1_density", passwd_)->set_initialized();
  }

  if (!S_->GetField("phase1_viscosity", passwd_)->initialized()) {
    *(S_->GetScalarData("phase1_viscosity", passwd_)) = 1.0;
    S_->GetField("phase1_viscosity", passwd_)->set_initialized();
  }

  if (!S_->GetField("phase2_density", passwd_)->initialized()) {
    *(S_->GetScalarData("phase2_density", passwd_)) = 1.0;
    S_->GetField("phase2_density", passwd_)->set_initialized();
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

  if (!S_->GetField(diff_tensor_name1_, passwd_)->initialized()) {
    S_->GetFieldData(diff_tensor_name1_, passwd_)->PutScalar(1.0);
    S_->GetField(diff_tensor_name1_, passwd_)->set_initialized();
  }

  if (!S_->GetField(diff_tensor_name2_, passwd_)->initialized()) {
    S_->GetFieldData(diff_tensor_name2_, passwd_)->PutScalar(1.0);
    S_->GetField(diff_tensor_name2_, passwd_)->set_initialized();
  }

  if (S_->HasField(molar_fraction_name1_)) {
    if (!S_->GetField(molar_fraction_name1_, passwd_)->initialized()) {
      S_->GetFieldData(molar_fraction_name1_, passwd_)->PutScalar(1.0);
      S_->GetField(molar_fraction_name1_, passwd_)->set_initialized();
    }
  }

  if (S_->HasField(molar_fraction_name2_)) {
    if (!S_->GetField(molar_fraction_name2_, passwd_)->initialized()) {
      S_->GetFieldData(molar_fraction_name2_, passwd_)->PutScalar(1.0);
      S_->GetField(molar_fraction_name2_, passwd_)->set_initialized();
    }
  }

  if (S_->HasField("fugacity1")) {
    if (!S_->GetField("fugacity1", passwd_)->initialized()) {
      S_->GetFieldData("fugacity1", passwd_)->PutScalar(1.0);
      S_->GetField("fugacity1", passwd_)->set_initialized();
    }
  }

  if (S_->HasField("fugacity2")) {
    if (!S_->GetField("fugacity2", passwd_)->initialized()) {
      S_->GetFieldData("fugacity2", passwd_)->PutScalar(1.0);
      S_->GetField("fugacity2", passwd_)->set_initialized();
    }
  }

  if (S_->HasField("saturation_w")) {
    if (!S_->GetField("saturation_w", passwd_)->initialized()) {
      S_->GetFieldData("saturation_w", passwd_)->PutScalar(0.0);
      S_->GetField("saturation_w", passwd_)->set_initialized();
    }
  }

  if (S_->HasField("saturation_n")) {
    if (!S_->GetField("saturation_n", passwd_)->initialized()) {
      S_->GetFieldData("saturation_n", passwd_)->PutScalar(0.0);
      S_->GetField("saturation_n", passwd_)->set_initialized();
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

  if (!S_->GetField("velocity_wet", passwd_)->initialized()) {
    S_->GetFieldData("velocity_wet", passwd_)->PutScalar(0.0);
    S_->GetField("velocity_wet", passwd_)->set_initialized();
  }

  if (!S_->GetField("velocity_nonwet", passwd_)->initialized()) {
    S_->GetFieldData("velocity_nonwet", passwd_)->PutScalar(0.0);
    S_->GetField("velocity_nonwet", passwd_)->set_initialized();
  } 
}


void Comp_PK::InitializeComponent()
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
  rho_ = *S_->GetScalarData("phase1_density");
  mu_ = *S_->GetScalarData("phase1_viscosity");
  phi_ = *S_->GetScalarData("porosity");

  // Initialize miscalleneous default parameters.
  ti_specs_ = NULL;
  error_control_ = FLOW_TI_ERROR_CONTROL_PRESSURE;

  src_sink_ = NULL;
  src_sink_distribution_ = 0;

  T_physics_ = 0.0;

  // create verbosity object
  vo_ = new VerboseObject("MPMC:Component", *comp_list_); 

  // Allocate memory for boundary data.
  bc_model_p_.resize(nfaces_wghost_, 0);
  bc_submodel_p_.resize(nfaces_wghost_, 0);
  bc_value_p_.resize(nfaces_wghost_, 0.0);
  bc_coef_p_.resize(nfaces_wghost_, 0.0);

  bc_model_s_.resize(nfaces_wghost_, 0);
  bc_submodel_s_.resize(nfaces_wghost_, 0);
  bc_value_s_.resize(nfaces_wghost_, 0.0);
  bc_coef_s_.resize(nfaces_wghost_, 0.0);

  bc_model_pc_.resize(nfaces_wghost_, 0);
  bc_submodel_pc_.resize(nfaces_wghost_, 0);
  bc_value_pc_.resize(nfaces_wghost_, 0.0);
  bc_coef_pc_.resize(nfaces_wghost_, 0.0);

  /*
  bc_model_cfw_.resize(nfaces_wghost_, 0);
  bc_submodel_cfw_.resize(nfaces_wghost_, 0);
  bc_value_cfw_.resize(nfaces_wghost_, 0.0);
  bc_coef_cfw_.resize(nfaces_wghost_, 0.0);

  bc_model_cfn_.resize(nfaces_wghost_, 0);
  bc_submodel_cfn_.resize(nfaces_wghost_, 0);
  bc_value_cfn_.resize(nfaces_wghost_, 0.0);
  bc_coef_cfn_.resize(nfaces_wghost_, 0.0);
  */

  bc_model_xw_.resize(nfaces_wghost_, 0);
  bc_submodel_xw_.resize(nfaces_wghost_, 0);
  bc_value_xw_.resize(nfaces_wghost_, 0.0);
  bc_coef_xw_.resize(nfaces_wghost_, 0.0);

  bc_model_xn_.resize(nfaces_wghost_, 0);
  bc_submodel_xn_.resize(nfaces_wghost_, 0);
  bc_value_xn_.resize(nfaces_wghost_, 0.0);
  bc_coef_xn_.resize(nfaces_wghost_, 0.0);
  bc_mixed_.resize(nfaces_wghost_, 0.0); 

  ProcessParameterList(*comp_list_);

  double time = S_->time();
  if (time >= 0.0) T_physics_ = time;

  time = T_physics_;
  bc_pressure_->Compute(time);
  //bc_flux_phase1_->Compute(time);
  //bc_flux_phase2_->Compute(time);
  bc_saturation_->Compute(time);
  bc_mfw_->Compute(time);
  bc_mfn_->Compute(time);

  // Create water retention models.
  rel_perm_w_ = Teuchos::rcp(new RelativePermeability(mesh_));
  rel_perm_n_ = Teuchos::rcp(new RelativePermeability(mesh_));
  rel_perm_w_->Init("wetting", wrm_list_);
  rel_perm_n_->Init("non wetting", wrm_list_);
  capillary_pressure_ = Teuchos::rcp(new CapillaryPressure(mesh_));
  capillary_pressure_->Init(wrm_list_); 
  capillary_pressure_old_ = Teuchos::rcp(new CapillaryPressure(mesh_));
  capillary_pressure_old_->Init(wrm_list_); 

  std::string krel_method_name = comp_list_->get<std::string>("relative permeability");
  //rel_perm_w_->ProcessStringRelativePermeability(krel_method_name);
  //rel_perm_n_->ProcessStringRelativePermeability(krel_method_name);

  ComputeBCs();
  ComputeBC_Pc();
  //ComputeBC_cf();
   
   /*
  //std::cout << "type of bc for mass fraction: \n";
  for (std::vector<int>::iterator i = bc_model_x_.begin(); i != bc_model_x_.end(); i++) {
    //std::cout << *i << "\n";
    if (*i == Operators::OPERATOR_BC_DIRICHLET) {
      bc_value_x_[i-bc_model_x_.begin()] = 0.5;
    }
  }
  */
  /*
  std::cout << "Value for bc of mass fraction: \n";
  for (std::vector<double>::iterator i = bc_value_x_.begin(); i != bc_value_x_.end(); i++) {
    std::cout << *i << "\n";
  }
  */

  // create operator boundary condition object
  op_bc_p_ = Teuchos::rcp(new Operators::BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model_p_, bc_value_p_, bc_mixed_));
  op_bc_pc_ = Teuchos::rcp(new Operators::BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model_pc_, bc_value_pc_, bc_mixed_));
  op_bc_xw_ = Teuchos::rcp(new Operators::BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model_xw_, bc_value_xw_, bc_mixed_));
  op_bc_xn_ = Teuchos::rcp(new Operators::BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model_xn_, bc_value_xn_, bc_mixed_));

  // allocate memory for absolute permeability
  K_.resize(ncells_wghost_);
  D1_.resize(ncells_wghost_);
  D2_.resize(ncells_wghost_);
  SetAbsolutePermeabilityTensor();
  SetDiffusionTensor();

  /*
  if (src_sink_distribution_ & Amanzi::CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
    Kxy_ = Teuchos::rcp(new Epetra_Vector(mesh_->cell_map(false)));
  }
  */

  // Select a proper matrix class. 
  Teuchos::ParameterList& tmp_list = comp_list_->sublist("operators").sublist("diffusion operator");
  Teuchos::ParameterList oplist_matrix = tmp_list.sublist("matrix");
  Teuchos::ParameterList oplist_pc = tmp_list.sublist("preconditioner");
  op_list_ = Teuchos::rcp(new Teuchos::ParameterList(tmp_list.sublist("preconditioner")));
  //std::cout << "op_list_: " << *op_list_ << "\n";

  // create diffusion operators op_matrix_ and op_preconditioner_.
  // They will need to be initialized, which is done later in InitNextTI.
  Operators::PDE_DiffusionFactory opfactory;
  op1_matrix_ = opfactory.Create(oplist_matrix, mesh_, op_bc_p_, 1.0, gravity_);
  op2_matrix_ = opfactory.Create(oplist_matrix, mesh_, op_bc_pc_, 1.0, gravity_);
  op3_matrix_ = opfactory.Create(oplist_matrix, mesh_, op_bc_xw_, 1.0, gravity_);
  op4_matrix_ = opfactory.Create(oplist_matrix, mesh_, op_bc_xn_, 1.0, gravity_);
  op1_preconditioner_ = opfactory.Create(oplist_pc, mesh_, op_bc_p_, 1.0, gravity_);
  op3_preconditioner_ = opfactory.Create(oplist_pc, mesh_, op_bc_xw_, 1.0, gravity_);

  for (int n_op = 0; n_op < soln_->size(); n_op++) {
    ops_.push_back(opfactory.Create(oplist_pc, mesh_, op_bc_p_, 1.0, gravity_));
  }

  Teuchos::ParameterList olist_adv = comp_list_->sublist("operators").sublist("advection operator");
  op2_preconditioner_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, mesh_));

  // initialize saturation of phase1 to compute flux of phase1 later
  //saturation_phase1_ = Teuchos::rcp(new CompositeVector(soln_->SubVector(1)->Data()->Map()));

  // init upwind operator
  Teuchos::ParameterList upw_list = comp_list_->sublist("operators")
                                            .sublist("diffusion operator")
                                            .sublist("upwind");

  //Operators::UpwindFactory<RelativePermeability> upwind_factory;
  //upwind_n_ = upwind_factory.Create(mesh_, rel_perm_n_, upw_list);
  //upwind_n_ = upwind_factory.Create(mesh_, rel_perm_w_, upw_list);

  upwind_n_ = Teuchos::rcp(new Operators::UpwindFlux<RelativePermeability>(mesh_, rel_perm_n_));
  upwind_n_->Init(upw_list);
  upwind_w_ = Teuchos::rcp(new Operators::UpwindFlux<RelativePermeability>(mesh_, rel_perm_w_));
  upwind_w_->Init(upw_list);
  //std::cout << "Done creating upwind operator \n";

  // init upwind_velocity
  upwind_vw_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("velocity_wet", passwd_)));
  upwind_vn_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("velocity_nonwet", passwd_)));
  //std::cout << "Done init upwind_velocity \n";

  tmp_flux1_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("velocity_wet", passwd_)));
  tmp_flux2_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("velocity_nonwet", passwd_)));
  //std::cout << "Done InitializeComponent\n";
} // End InitializePressure()


/* ******************************************************************
* Generic initialization of a next time integration phase.
****************************************************************** */
void Comp_PK::InitNextTI()
{
  // initialize operators
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > D1ptr = Teuchos::rcpFromRef(D1_);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > D2ptr = Teuchos::rcpFromRef(D2_);
  op1_matrix_->Setup(Kptr, rel_perm_w_->Krel(), Teuchos::null);
  op1_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op2_matrix_->Setup(Kptr, rel_perm_n_->Krel(), Teuchos::null);
  op2_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op3_matrix_->Setup(D1ptr, Teuchos::null, Teuchos::null);
  op3_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op4_matrix_->Setup(D2ptr, Teuchos::null, Teuchos::null);
  op4_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op1_preconditioner_->Setup(Kptr, rel_perm_w_->Krel(), Teuchos::null);
  op1_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op2_preconditioner_->Setup(*upwind_vw_);
  op2_preconditioner_->UpdateMatrices(*upwind_vw_);

  for (int ii = 0; ii < ops_.size(); ii++) {
    ops_[ii]->Setup(Kptr, rel_perm_n_->Krel(), Teuchos::null);
    ops_[ii]->UpdateMatrices(Teuchos::null, Teuchos::null);
  }
} // End InitNextTI


/* ******************************************************************* 
* Wrapper for Comp_PK::CommitState. This is the interface used 
* by the MPC.
******************************************************************* */
void Comp_PK::CommitStep(double t_old, double t_new)
{
  // Write pressure to state
  CompositeVector& p2 = *S_->GetFieldData("pressure_n", passwd_);
  CompositeVector& p1 = *S_->GetFieldData("pressure_w", passwd_);
  p1 = *soln_->SubVector(0)->Data();

  // write saturation of phase1 to state
  CompositeVector& s2 = *S_->GetFieldData("saturation_n", passwd_);
  s2 = *soln_->SubVector(1)->Data();
  CompositeVector& s1 = *S_->GetFieldData("saturation_w", passwd_);
  s1.Update(-1.0, s2, 0.0);
  s1.Shift(1.0);

  capillary_pressure_->Compute(s1); 
  p2.Update(1.0, p1, 0.0);
  p2.Update(1.0, *capillary_pressure_->Pc(), 1.0);

  // compute upwind velocities from pressure 
  CompositeVector& phase1_flux = *S_->GetFieldData("velocity_wet", passwd_);
  op1_matrix_->UpdateFlux(p1, phase1_flux);
  CompositeVector& phase2_flux = *S_->GetFieldData("velocity_nonwet", passwd_);
  op2_matrix_->UpdateFlux(p2, phase2_flux);

  // compute molar fractions
  const CompositeVector& fuga_comp = *soln_->SubVector(1+comp_id_)->Data();
  CompositeVector& xw = *S_->GetFieldData(molar_fraction_name1_, passwd_);
  CompositeVector& xn = *S_->GetFieldData(molar_fraction_name2_, passwd_);
  xn.ReciprocalMultiply(1.0, p2, fuga_comp, 0.0);
  xw.Update(1.0/comp_coeff_, fuga_comp, 0.0);
}


/* ******************************************************************
*  Temporary convertion from double to tensor.                                               
****************************************************************** */
void Comp_PK::SetAbsolutePermeabilityTensor()
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


/* ******************************************************************
*  Temporary convertion from double to tensor.                                               
****************************************************************** */
void Comp_PK::SetDiffusionTensor()
{
  const CompositeVector& cv1 = *S_->GetFieldData(diff_tensor_name1_);
  cv1.ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& perm1 = *cv1.ViewComponent("cell", true);

  if (dim_ == 2) {
    for (int c = 0; c < D1_.size(); c++) {
      if (perm1[0][c] == perm1[1][c]) {
        D1_[c].Init(dim_, 1);
        D1_[c](0, 0) = perm1[0][c];
      } else {
        D1_[c].Init(dim_, 2);
        D1_[c](0, 0) = perm1[0][c];
        D1_[c](1, 1) = perm1[1][c];
      }
    }    
  } else if (dim_ == 3) {
    for (int c = 0; c < D1_.size(); c++) {
      if (perm1[0][c] == perm1[1][c] && perm1[0][c] == perm1[2][c]) {
        D1_[c].Init(dim_, 1);
        D1_[c](0, 0) = perm1[0][c];
      } else {
        D1_[c].Init(dim_, 2);
        D1_[c](0, 0) = perm1[0][c];
        D1_[c](1, 1) = perm1[1][c];
        D1_[c](2, 2) = perm1[2][c];
      }
    }        
  }

  const CompositeVector& cv2 = *S_->GetFieldData(diff_tensor_name2_);
  cv2.ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& perm2 = *cv2.ViewComponent("cell", true);

  if (dim_ == 2) {
    for (int c = 0; c < D2_.size(); c++) {
      if (perm2[0][c] == perm2[1][c]) {
        D2_[c].Init(dim_, 1);
        D2_[c](0, 0) = perm2[0][c];
      } else {
        D2_[c].Init(dim_, 2);
        D2_[c](0, 0) = perm2[0][c];
        D2_[c](1, 1) = perm2[1][c];
      }
    }    
  } else if (dim_ == 3) {
    for (int c = 0; c < D2_.size(); c++) {
      if (perm2[0][c] == perm2[1][c] && perm2[0][c] == perm2[2][c]) {
        D2_[c].Init(dim_, 1);
        D2_[c](0, 0) = perm2[0][c];
      } else {
        D2_[c].Init(dim_, 2);
        D2_[c](0, 0) = perm2[0][c];
        D2_[c](1, 1) = perm2[1][c];
        D2_[c](2, 2) = perm2[2][c];
      }
    }        
  }
}


// need to fix flag essential bc. We need Dirichlet bc for all 4 primary
// variables before incrementing Dirichlet bc counter 
/* ******************************************************************* 
* Compute boundary conditions based on bc type and values
******************************************************************* */
void Comp_PK::ComputeBCs() {
  int flag_essential_bc = 0;
  dirichlet_bc_faces_ = 0;

  for (int n = 0; n < bc_model_p_.size(); n++) {
    bc_model_p_[n] = Operators::OPERATOR_BC_NONE;
    bc_value_p_[n] = 0.0;
    bc_coef_p_[n] = 0.0;
    bc_model_s_[n] = Operators::OPERATOR_BC_NONE;
    bc_value_s_[n] = 0.0;
    bc_coef_s_[n] = 0.0;
    bc_model_xw_[n] = Operators::OPERATOR_BC_NONE;
    bc_value_xw_[n] = 0.0;
    bc_coef_xw_[n] = 0.0;   
    bc_model_xn_[n] = Operators::OPERATOR_BC_NONE;
    bc_value_xn_[n] = 0.0;
    bc_coef_xn_[n] = 0.0;   
  }

  Flow::FlowBoundaryFunction::Iterator bc;
  for (bc = bc_pressure_->begin(); bc != bc_pressure_->end(); ++bc) {
    int f = bc->first;
    bc_model_p_[f] = Operators::OPERATOR_BC_DIRICHLET;
    bc_value_p_[f] = bc->second;
    flag_essential_bc = 1;
    dirichlet_bc_faces_++;
  }

  for (bc = bc_saturation_->begin(); bc != bc_saturation_->end(); ++bc) {
    int f = bc->first;
    bc_model_s_[f] = Operators::OPERATOR_BC_DIRICHLET;
    bc_value_s_[f] = 1.0 - bc->second;
  }

  for (bc = bc_mfw_->begin(); bc != bc_mfw_->end(); ++bc) {
    int f = bc->first;
    bc_model_xw_[f] = Operators::OPERATOR_BC_DIRICHLET;
    bc_value_xw_[f] = bc->second;
  }

  for (bc = bc_mfn_->begin(); bc != bc_mfn_->end(); ++bc) {
    int f = bc->first;
    bc_model_xn_[f] = Operators::OPERATOR_BC_DIRICHLET;
    bc_value_xn_[f] = bc->second;
  }

  /*
  for (bc = bc_flux_phase1_->begin(); bc != bc_flux_phase1_->end(); ++bc) {
    int f = bc->first;
    bc_model_p_[f] = Operators::OPERATOR_BC_NEUMANN;
    bc_value_p_[f] = bc->second;
  }

  for (bc = bc_flux_phase2_->begin(); bc != bc_flux_phase2_->end(); ++bc) {
    int f = bc->first;
    bc_model_s_[f] = Operators::OPERATOR_BC_NEUMANN;
    bc_value_s_[f] = bc->second;
  }
  */

  // mark missing boundary conditions as zero flux conditions 
  AmanziMesh::Entity_ID_List cells;
  missed_bc_faces_ = 0;
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_model_p_[f] == Operators::OPERATOR_BC_NONE) {
      cells.clear();
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int ncells = cells.size();

      if (ncells == 1) {
        bc_model_p_[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value_p_[f] = 0.0;
        bc_model_s_[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value_s_[f] = 0.0;
        bc_model_xw_[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value_xw_[f] = 0.0;
        bc_model_xn_[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value_xn_[f] = 0.0;
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
void Comp_PK::ComputeBC_Pc()
{
  //std::cout << "Compute boundary condition for Pc\n";
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm_n_->WRM(); 
  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    //std::cout << "block " << mb << " region: " << region << "\n";
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED, &block);
    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      //std::cout << "Processing boundary face " << *i << "\n";
      if (bc_model_p_[*i] == Operators::OPERATOR_BC_DIRICHLET) {
        bc_model_pc_[*i] = Operators::OPERATOR_BC_DIRICHLET;
        //std::cout << "Dirichlet value of saturation on face " << *i << " = " << bc_value_s_[*i] << "\n";
        bc_value_pc_[*i] = bc_value_p_[*i] + WRM[mb]->capillaryPressure(bc_value_s_[*i]);
        //std::cout << "phase2 pressure on face " << *i << " = " << bc_value_pc_[*i] << "\n";
      } else {
        bc_model_pc_[*i] = bc_model_p_[*i];
        bc_value_pc_[*i] = bc_value_p_[*i];
        // do nothing for now
      }
    }
  }
  //std::cout << "done computing capillary pressure \n";
}

void Comp_PK::DeriveFaceValuesFromCellValues(
    const Epetra_MultiVector& ucells, Epetra_MultiVector& ufaces,
    const std::vector<int>& bc_model, const std::vector<double>& bc_value)
{
  AmanziMesh::Entity_ID_List cells;
  int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  for (int f = 0; f < nfaces; f++) {
    cells.clear();
    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::OWNED, &cells);
    int ncells = cells.size();

    if (ncells == 2) {
      double face_value = 0.0;
      for (int n = 0; n < ncells; n++) face_value += ucells[0][cells[n]];
      ufaces[0][f] = face_value / ncells;
    } else if (ncells == 1) {
        if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
          ufaces[0][f] = bc_value[f];
        } else {
          ufaces[0][f] = ucells[0][cells[0]];
        }
    } else {
      // do nothing
    }
  }
}


/* ******************************************************************
* Compute boundary conditions for capillary pressure                                          
****************************************************************** */
/*
void Comp_PK::ComputeBC_cf()
{
  //std::cout << "Compute boundary condition for Pc\n";
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm_n_->WRM(); 
  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    //std::cout << "block " << mb << " region: " << region << "\n";
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED, &block);
    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      //std::cout << "Processing boundary face " << *i << "\n";
      if (bc_model_s_[*i] == Operators::OPERATOR_BC_DIRICHLET) {
        bc_model_cfw_[*i] = bc_model_s_[*i];
        bc_model_cfn_[*i] = bc_model_s_[*i];
        //std::cout << "Dirichlet value of saturation on face " << *i << " = " << bc_value_s_[*i] << "\n";
        bc_value_cfw_[*i] = bc_value_xw_[*i] * WRM[mb]->k_relative(1.0 - bc_value_s_[*i], "wetting");
        bc_value_cfn_[*i] = bc_value_xn_[*i] * WRM[mb]->k_relative(1.0 - bc_value_s_[*i], "non wetting");
        //std::cout << "phase2 pressure on face " << *i << " = " << bc_value_pc_[*i] << "\n";
      } else {
        bc_model_cfw_[*i] = bc_model_s_[*i];
        bc_value_cfw_[*i] = bc_value_s_[*i];
        bc_model_cfn_[*i] = bc_model_s_[*i];
        bc_value_cfn_[*i] = bc_value_s_[*i];
        // do nothing for now
      }
    }
  }
  //std::cout << "done computing capillary pressure \n";
}
*/
}  // namespace Flow
} // End namespace Amanzi
