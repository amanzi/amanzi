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
#include "OperatorDefs.hh"

#include "Mesh.hh"

#include "Phase2_PK.hh"
#include "MultiphaseDefs.hh"

// namespace
namespace Amanzi {
namespace Multiphase {

//Simple construction for multiphase pressure pk
Phase2_PK::Phase2_PK(Teuchos::ParameterList& pk_tree,
                      const Teuchos::RCP<Teuchos::ParameterList>& glist,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& soln):
              bc_pressure_(NULL),
              bc_flux_phase2_(NULL),
              bc_saturation_(NULL),
              src_sink_(NULL),
              S_(S),
              soln_(soln),
              passwd_("state"),
              vo_(NULL),
              include_capillary_(false)
{
  glist_ = Teuchos::rcp(new Teuchos::ParameterList(*glist));
  mesh_ = S->GetMesh();
  dim_ = mesh_->space_dimension();
  src_sink_ = NULL;
  src_sink_distribution_ = 0;
    
  // We need the flow list
  Teuchos::ParameterList flow_list;
  if (glist->isSublist("Flow")) {
    flow_list = glist->sublist("Flow");
  } else {
    Errors::Message msg("Multiphase Phase2 PK: input parameter list does not have <Flow> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  if (flow_list.isSublist("Phase2 Specs")) {
    mp_list_ = flow_list.sublist("Phase2 Specs");
  } else {
    Errors::Message msg("Multiphase Phase2 PK: input parameter list does not have <Phase2 Specs> sublist.");
    Exceptions::amanzi_throw(msg);
  }

  // require state variables for the Multiphase Phase2
  if (!S->HasField("phase2_density")) {
    S->RequireScalar("phase2_density", passwd_);
  }
  if (!S->HasField("phase2_viscosity")) {
    S->RequireScalar("phase2_viscosity", passwd_);
  }

  if (!S->HasField("gravity")) {
    S->RequireConstantVector("gravity", passwd_, dim_);
  }

  if (!S->HasField("porosity")) {
    S->RequireField("porosity", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField("phase2_pressure")) {
    S->RequireField("phase2_pressure", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  /*
  if (!S->HasField("perm_x")) {
    S->RequireField("perm_x", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S->HasField("perm_y")) {
    S->RequireField("perm_y", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S->HasField("perm_z")) {
    S->RequireField("perm_z", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  */
  
  if (!S->HasField("permeability")) {
    S->RequireField("permeability", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim_);
  }

  if (!S->HasField("phase2_saturation")) {
    S->RequireField("phase2_saturation", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField("phase2_velocity")) {
    S->RequireField("phase2_velocity", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);
  }
} // End Phase2_PK constructor


/* ******************************************************************
* Clean memory.
****************************************************************** */
Phase2_PK::~Phase2_PK()
{
  if (bc_pressure_ != NULL) delete bc_pressure_;
  if (bc_flux_phase2_ != NULL) delete bc_flux_phase2_;

  if (src_sink_ != NULL) delete src_sink_;
  if (vo_ != NULL) delete vo_;
}


void Phase2_PK::InitializeFields()
{
  // set popular default values

  if (!S_->GetField("phase2_density", passwd_)->initialized()) {
    *(S_->GetScalarData("phase2_density", passwd_)) = 1.0;
    S_->GetField("phase2_density", passwd_)->set_initialized();
  }

  if (!S_->GetField("phase2_viscosity", passwd_)->initialized()) {
    *(S_->GetScalarData("phase2_viscosity", passwd_)) = 1.0;
    S_->GetField("phase2_viscosity", passwd_)->set_initialized();
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

  /*
  if (!S_->GetField("perm_x", passwd_)->initialized()) {
    S_->GetFieldData("perm_x", passwd_)->PutScalar(1.0);
    S_->GetField("perm_x", passwd_)->set_initialized();
  }
  if (!S_->GetField("perm_y", passwd_)->initialized()) {
    S_->GetFieldData("perm_y", passwd_)->PutScalar(1.0);
    S_->GetField("perm_y", passwd_)->set_initialized();
  }
  if (!S_->GetField("perm_z", passwd_)->initialized()) {
    S_->GetFieldData("perm_z", passwd_)->PutScalar(1.0);
    S_->GetField("perm_z", passwd_)->set_initialized();
  }   
  */

  if (!S_->GetField("porosity", passwd_)->initialized()) {
    S_->GetFieldData("porosity", passwd_)->PutScalar(1.0);
    S_->GetField("porosity", passwd_)->set_initialized();
  }

  if (!S_->GetField("phase2_saturation", passwd_)->initialized()) {
    S_->GetFieldData("phase2_saturation", passwd_)->PutScalar(0.0);
    S_->GetField("phase2_saturation", passwd_)->set_initialized();
  }

  if (!S_->GetField("phase2_pressure", passwd_)->initialized()) {
    S_->GetFieldData("phase2_pressure", passwd_)->PutScalar(0.0);
    S_->GetField("phase2_pressure", passwd_)->set_initialized();
  }

  if (!S_->GetField("phase2_velocity", passwd_)->initialized()) {
    S_->GetFieldData("phase2_velocity", passwd_)->PutScalar(0.0);
    S_->GetField("phase2_velocity", passwd_)->set_initialized();
  }
}


void Phase2_PK::InitializePhase2()
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
  rho_ = *S_->GetScalarData("phase2_density");
  mu_ = *S_->GetScalarData("phase2_viscosity");

  //rho_ = 1.0;

  src_sink_ = NULL;
  src_sink_distribution_ = 0;

  // create verbosity object
  vo_ = new VerboseObject("MultiphasePK::Phase2", mp_list_); 

  ProcessParameterList(mp_list_);

  // compute boundary values from input functions
  T_physics_ = 0.0;
  double time = S_->time();
  if (time >= 0.0) T_physics_ = time;

  time = T_physics_;
  bc_pressure_->Compute(time, time);
  bc_flux_phase2_->Compute(time, time);
  bc_saturation_->Compute(time, time);

  // allocate memory for absolute permeability
  K_.resize(ncells_wghost_);

  if (src_sink_distribution_ & Amanzi::CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
    Kxy_ = Teuchos::rcp(new Epetra_Vector(mesh_->cell_map(false)));
  }
  SetAbsolutePermeabilityTensor();
  
  // Create water retention models.
  if (!mp_list_.isSublist("Water retention models")) {
    Errors::Message msg;
    msg << "Multiphase PK: input XML files does not have <Water retention models> sublist in <Multiphase problem> \n";
    Exceptions::amanzi_throw(msg);
  }
  Teuchos::ParameterList& wrm_list = mp_list_.sublist("Water retention models");
  rel_perm_n_ = Teuchos::rcp(new RelativePermeability(mesh_));
  rel_perm_n_->Init("non wetting", wrm_list);
  capillary_pressure_ = Teuchos::rcp(new CapillaryPressure(mesh_));
  capillary_pressure_->Init(wrm_list); 

  std::string krel_method_name = mp_list_.get<std::string>("relative permeability");
  // rel_perm_n_->ProcessStringRelativePermeability(krel_method_name);

  // Compute boundary conditions
  ComputeBCs();
  ComputeBC_Pc();
  
  // create operator boundary condition object
  op_bc_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  op_bc_s_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  op_bc_pc_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  op_bc_pc_prime_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

  // Select a proper matrix class. 
  Teuchos::ParameterList& tmp_list = mp_list_.sublist("operators").sublist("diffusion operator");
  Teuchos::ParameterList oplist_matrix = tmp_list.sublist("matrix");
  Teuchos::ParameterList oplist_pc = tmp_list.sublist("preconditioner");

  // create diffusion operators op_matrix_ and op_preconditioner_.
  // They will need to be initialized, which is done later in InitNextTI.
  Operators::PDE_DiffusionFactory opfactory;
  op_matrix_ = opfactory.Create(oplist_matrix, mesh_, op_bc_pc_, rho_, gravity_);
  //op_matrix_->SetDensity(rho_);
  op_matrix_copy_ = opfactory.Create(oplist_matrix, mesh_, op_bc_pc_prime_, rho_, gravity_);
  //op_matrix_copy_->SetDensity(rho_);
  op_matrix_copy1_ = opfactory.Create(oplist_matrix, mesh_, op_bc_pc_, rho_, gravity_);
  //op_matrix_copy1_->SetDensity(rho_);
  op1_preconditioner_ = opfactory.Create(oplist_pc, mesh_, op_bc_, rho_, gravity_);
  for (int n_op = 0; n_op < soln_->size(); n_op++) {
    ops_.push_back(opfactory.Create(oplist_pc, mesh_, op_bc_pc_, rho_, gravity_));
  }

  Teuchos::ParameterList olist_adv = mp_list_.sublist("operators").sublist("advection operator");
  op2_preconditioner_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, mesh_));
  op_sum1_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, op2_preconditioner_->global_operator()));
  op_sum_ = Teuchos::rcp(new Operators::PDE_DiffusionFV(oplist_pc, op_sum1_->global_operator()));

  // accumulation operators
  op_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op2_preconditioner_->global_operator())); 

  // create pressure_phase2_ which is needed for Functional and UpdatePreconditioner
  pressure_phase2_ = Teuchos::rcp(new CompositeVector(soln_->SubVector(0)->Data()->Map()));

  // init upwind operator
  Teuchos::ParameterList upw_list = mp_list_.sublist("operators")
                                            .sublist("diffusion operator")
                                            .sublist("upwind");
  Operators::UpwindFactory<RelativePermeability> upwind_factory;
  upwind_ = upwind_factory.Create(mesh_, rel_perm_n_, upw_list);
  Operators::UpwindFactory<CapillaryPressure> upwind_pc_factory;
  upwind_pc_ = upwind_pc_factory.Create(mesh_, capillary_pressure_, upw_list);

  // init upwind_velocity
  upwind_velocity_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("phase2_velocity", passwd_)));
  tmp_flux_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("phase2_velocity", passwd_)));

  //std::cout << "Done InitializePhase2\n";
} // End InitializePressure()


/* ******************************************************************
* Generic initialization of a next time integration phase.
****************************************************************** */
void Phase2_PK::InitNextTI()
{
  std::vector<int>& bc_model_s = op_bc_s_->bc_model();

  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);
  cvs.AddComponent("face", AmanziMesh::FACE, 1);

  Teuchos::RCP<CompositeVector> Krel = Teuchos::rcp(new CompositeVector(cvs));
  Krel->PutScalarMasterAndGhosted(0.0);

  // Calculate phase relative permeability needed to initialize diffusion operator
  CompositeVector& phase2_saturation = *S_->GetFieldData("phase2_saturation", passwd_);
  Teuchos::RCP<CompositeVector> phase1_sat = Teuchos::rcp(new CompositeVector(phase2_saturation));
  phase1_sat->Scale(-1.0);
  phase1_sat->Shift(1.0);
  rel_perm_n_->Compute(*phase1_sat);
  rel_perm_n_->Krel()->Scale(rho_/mu_);
  rel_perm_n_->dKdS()->Scale(-rho_/mu_);

  upwind_velocity_ = S_->GetFieldData("phase2_velocity", passwd_);
  upwind_->Compute(*upwind_velocity_, *upwind_velocity_, bc_model_s_, *rel_perm_n_->Krel());
  upwind_->Compute(*upwind_velocity_, *upwind_velocity_, bc_model_s_, *rel_perm_n_->dKdS());

  // initialize matrix and preconditioner operators 
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_); 
  op_matrix_->Setup(Kptr, rel_perm_n_->Krel(), Teuchos::null);
  op_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op1_preconditioner_->Setup(Kptr, rel_perm_n_->Krel(), Teuchos::null);
  op1_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // this operator needs newton correction
  op2_preconditioner_->Setup(*upwind_velocity_);
  op2_preconditioner_->UpdateMatrices(rel_perm_n_->dKdS().ptr(), Teuchos::null);

  for (int ii = 0; ii < ops_.size(); ii++) {
    ops_[ii]->Setup(Kptr, rel_perm_n_->Krel(), Teuchos::null);
    ops_[ii]->UpdateMatrices(Teuchos::null, Teuchos::null);
  }

  /*
  // Initialize source
  if (src_sink_ != NULL) {
    if (src_sink_distribution_ & Amanzi::CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      src_sink_->ComputeDistribute(T0, Kxy_->Values()); 
    } else {
      src_sink_->ComputeDistribute(T0, NULL);
    }
  }
  */
} // End InitNextTI


bool Phase2_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  // do nothing
}


/* ******************************************************************* 
* Wrapper for Phase1_PK::CommitState. This is the interface used 
* by the MPC.
******************************************************************* */
void Phase2_PK::CommitStep(double t_old, double t_new)
{
  // Write pressure of phase1 to state
  CompositeVector& p2 = *S_->GetFieldData("phase2_pressure", passwd_);
  const CompositeVector& p1 = *soln_->SubVector(0)->Data();

  // write saturation of phase1 to state
  CompositeVector& s2 = *S_->GetFieldData("phase2_saturation", passwd_);
  s2 = *soln_->SubVector(1)->Data();

  // compute P2
  if (include_capillary_) {
    Teuchos::RCP<CompositeVector> s1 = Teuchos::rcp(new CompositeVector(s2));
    s1->Scale(-1.0);
    s1->Shift(1.0);
    capillary_pressure_->Compute(*s1); 
    p2.Update(1.0, p1, 0.0);
    p2.Update(1.0, *capillary_pressure_->Pc(), 1.0);
  } else {
    p2.Update(1.0, p1, 0.0); 
  }

  // compute phase2 flux from p2
  CompositeVector& phase2_flux = *S_->GetFieldData("phase2_velocity", passwd_);
  op_matrix_->UpdateFlux(p2, phase2_flux);
}


/* ******************************************************************
* Add source and sink terms.                                   
****************************************************************** */
void Phase2_PK::AddSourceTerms(CompositeVector& rhs)
{
  Epetra_MultiVector& rhs_cell = *rhs.ViewComponent("cell", true);
  PK_DomainFunction::Iterator src;
  for (src = src_sink_->begin(); src != src_sink_->end(); ++src) {
    int c = src->first;
    rhs_cell[0][c] += mesh_->cell_volume(c) * src->second[0];
  }
}


/* ******************************************************************
*  Temporary convertion from double to tensor.                                               
****************************************************************** */
void Phase2_PK::SetAbsolutePermeabilityTensor()
{
  const CompositeVector& cv = *S_->GetFieldData("permeability");
  cv.ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& perm = *cv.ViewComponent("cell", true);
  /*
  const CompositeVector& cvx = *S_->GetFieldData("perm_x");
  cvx.ScatterMasterToGhosted("cell");
  const CompositeVector& cvy = *S_->GetFieldData("perm_y");
  cvy.ScatterMasterToGhosted("cell");
  const CompositeVector& cvz = *S_->GetFieldData("perm_z");
  cvz.ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& permx = *cvx.ViewComponent("cell", true);
  const Epetra_MultiVector& permy = *cvy.ViewComponent("cell", true);
  const Epetra_MultiVector& permz = *cvz.ViewComponent("cell", true);
  */

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


/* ******************************************************************* 
* Compute boundary conditions based on bc type and values
******************************************************************* */
void Phase2_PK::ComputeBCs()
{
  std::vector<int>& bc_model_p = op_bc_->bc_model();

  int flag_essential_bc = 0;
  dirichlet_bc_faces_ = 0;

  for (int n = 0; n < bc_model_p_.size(); n++) {
    bc_model_p_[n] = Operators::OPERATOR_BC_NONE;
    bc_value_p_[n] = 0.0;
    bc_coef_p_[n] = 0.0;
    bc_model_s_[n] = Operators::OPERATOR_BC_NONE;
    bc_value_s_[n] = 0.0;
    bc_coef_s_[n] = 0.0;
  }

  Flow::FlowBoundaryFunction::Iterator bc;
  for (bc = bc_pressure_->begin(); bc != bc_pressure_->end(); ++bc) {
    int f = bc->first;
    bc_model_p_[f] = Operators::OPERATOR_BC_DIRICHLET;
    bc_value_p_[f] = bc->second[0];
    flag_essential_bc = 1;
    dirichlet_bc_faces_++;
  }

  for (bc = bc_saturation_->begin(); bc != bc_saturation_->end(); ++bc) {
    int f = bc->first;
    // currently hard-coded s2 as primary variable.
    // However, rel perm and capillary pressure are computed using s1
    // this only matters when we need to compute the derivatives of Krel and Pc
    // since dS2 = - dS1.
    bc_model_s_[f] = Operators::OPERATOR_BC_DIRICHLET;
    bc_value_s_[f] = 1.0 - bc->second[0]; // bc value for s1
  }

  for (bc = bc_flux_phase2_->begin(); bc != bc_flux_phase2_->end(); ++bc) {
    int f = bc->first;
    bc_model_p_[f] = Operators::OPERATOR_BC_NEUMANN;
    bc_value_p_[f] = bc->second[0];
  }
  
  // mark missing boundary conditions as zero flux conditions only for pressure equation
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
void Phase2_PK::ComputeBC_Pc()
{
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm_n_->WRM(); 
  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED, &block);
    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      //std::cout << "Processing boundary face " << *i << "\n";
      if (bc_model_p_[*i] == Operators::OPERATOR_BC_DIRICHLET) {
        bc_model_pc_[*i] = Operators::OPERATOR_BC_DIRICHLET;
        bc_model_pc_prime_[*i] = Operators::OPERATOR_BC_DIRICHLET;
        //std::cout << "Dirichlet value of saturation on face " << *i << " = " << bc_value_s_[*i] << "\n";
        bc_value_pc_[*i] = bc_value_p_[*i] + WRM[mb]->capillaryPressure(bc_value_s_[*i]); // compute Pc using s1
        //std::cout << "phase2 pressure on face " << *i << " = " << bc_value_pc_[*i] << "\n";
        bc_value_pc_prime_[*i] = -WRM[mb]->dPc_dS(bc_value_s_[*i]);
      } else {
        bc_model_pc_[*i] = bc_model_p_[*i];
        bc_value_pc_[*i] = bc_value_p_[*i];

        bc_model_pc_prime_[*i] = bc_model_p_[*i];
        bc_value_pc_prime_[*i] = bc_value_p_[*i];
      }
    }
  }
}

}  // namespace Flow
} // End namespace Amanzi
