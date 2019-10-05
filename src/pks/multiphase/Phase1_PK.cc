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

#include "Phase1_PK.hh"
#include "MultiphaseDefs.hh"

// namespace
namespace Amanzi {
namespace Multiphase {

//Simple construction for multiphase pressure pk
Phase1_PK::Phase1_PK(Teuchos::ParameterList& pk_tree,
                      const Teuchos::RCP<Teuchos::ParameterList>& glist,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& soln):
              bc_pressure_(NULL),
              bc_flux_phase1_(NULL),
              bc_saturation_(NULL),
              src_sink_(NULL),
              S_(S),
              soln_(soln),
              passwd_("state"),
              vo_(NULL)
{
  mesh_ = S->GetMesh();
  dim_ = mesh_->space_dimension();
  src_sink_ = NULL;
  src_sink_distribution_ = 0;
    
  // We need the flow list
  Teuchos::RCP<Teuchos::ParameterList> flow_list = Teuchos::sublist(glist, "Flow", true);
  mp_list_ = Teuchos::sublist(flow_list, "Phase1 Specs", true);

  // require state variables for the MultiphasePhase1 
  if (!S->HasField("phase1_density")) {
    S->RequireScalar("phase1_density", passwd_);
  }
  if (!S->HasField("phase1_viscosity")) {
    S->RequireScalar("phase1_viscosity", passwd_);
  }
  /*
  if (!S->HasField("porosity")) {
    S->RequireScalar("porosity", passwd_);
  }
  */

  if (!S->HasField("gravity")) {
    S->RequireConstantVector("gravity", passwd_, dim_);
  }

  if (!S->HasField("porosity")) {
    S->RequireField("porosity", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField("phase1_pressure")) {
    S->RequireField("phase1_pressure", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  /*
  if (!S->HasField("perm_x")) {
    S->RequireField("perm_x", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1.0);
  }
  if (!S->HasField("perm_y")) {
    S->RequireField("perm_y", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1.0);
  }
  if (!S->HasField("perm_z")) {
    S->RequireField("perm_z", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1.0);
  }
  */
  
  if (!S->HasField("permeability")) {
    S->RequireField("permeability", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim_);
  }

  if (!S->HasField("phase1_rel_perm")) {
    S->RequireField("phase1_rel_perm", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);
  }

  if (!S->HasField("phase1_saturation")) {
    S->RequireField("phase1_saturation", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField("phase1_velocity")) {
    S->RequireField("phase1_velocity", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("face", AmanziMesh::FACE, 1);
  }
} // End Phase1_PK constructor


/* ******************************************************************
* Clean memory.
****************************************************************** */
Phase1_PK::~Phase1_PK()
{
  if (bc_pressure_ != NULL) delete bc_pressure_;
  if (bc_flux_phase1_ != NULL) delete bc_flux_phase1_;

  if (src_sink_ != NULL) delete src_sink_;
  if (vo_ != NULL) delete vo_;
}


void Phase1_PK::InitializeFields()
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

  if (!S_->GetField("porosity", passwd_)->initialized()) {
    S_->GetFieldData("porosity", passwd_)->PutScalar(1.0);
    S_->GetField("porosity", passwd_)->set_initialized();
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

  if (!S_->GetField("phase1_saturation", passwd_)->initialized()) {
    S_->GetFieldData("phase1_saturation", passwd_)->PutScalar(0.0);
    S_->GetField("phase1_saturation", passwd_)->set_initialized();
  }

  if (!S_->GetField("phase1_pressure", passwd_)->initialized()) {
    S_->GetFieldData("phase1_pressure", passwd_)->PutScalar(0.0);
    S_->GetField("phase1_pressure", passwd_)->set_initialized();
  }

  if (!S_->GetField("phase1_rel_perm", passwd_)->initialized()) {
    S_->GetFieldData("phase1_rel_perm", passwd_)->PutScalar(0.0);
    S_->GetField("phase1_rel_perm", passwd_)->set_initialized();
  }

  if (!S_->GetField("phase1_velocity", passwd_)->initialized()) {
    S_->GetFieldData("phase1_velocity", passwd_)->PutScalar(0.0);
    S_->GetField("phase1_velocity", passwd_)->set_initialized();
  }
}


void Phase1_PK::InitializePhase1()
{
  // Initilize various common data depending on mesh and state.
  ncells_owned_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  ncells_wghost_ = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

  nfaces_owned_ = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nfaces_wghost_ = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

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
  //phi_ = *S_->GetScalarData("porosity");

  src_sink_ = NULL;
  src_sink_distribution_ = 0;

  // create verbosity object
  vo_ = new VerboseObject("MP::Phase1", *mp_list_); 

  // Allocate memory for boundary data.
  bc_model_p_.resize(nfaces_wghost_, 0);
  bc_submodel_p_.resize(nfaces_wghost_, 0);
  bc_value_p_.resize(nfaces_wghost_, 0.0);
  bc_coef_p_.resize(nfaces_wghost_, 0.0);

  bc_model_s_.resize(nfaces_wghost_, 0);
  bc_submodel_s_.resize(nfaces_wghost_, 0);
  bc_value_s_.resize(nfaces_wghost_, 0.0);
  bc_coef_s_.resize(nfaces_wghost_, 0.0);
  bc_mixed_.resize(nfaces_wghost_, 0.0); 

  ProcessParameterList(*mp_list_);

  bc_pressure_->Compute(0.0);
  bc_flux_phase1_->Compute(0.0);
  bc_saturation_->Compute(0.0);
  ComputeBCs();

  // create operator boundary condition object
  op_bc_p_ = Teuchos::rcp(new Operators::BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model_p_, bc_value_p_, bc_mixed_));
  op_bc_s_ = Teuchos::rcp(new Operators::BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_model_s_, bc_value_s_, bc_mixed_));

  // allocate memory for absolute permeability
  K_.resize(ncells_wghost_);

  if (src_sink_distribution_ & Amanzi::CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
    Kxy_ = Teuchos::rcp(new Epetra_Vector(mesh_->cell_map(false)));
  }
  SetAbsolutePermeabilityTensor();

  // Create water retention models.
  Teuchos::RCP<Teuchos::ParameterList> wrm_list = Teuchos::sublist(mp_list_, "Water retention models", true);
  rel_perm_w_ = Teuchos::rcp(new RelativePermeability(mesh_));
  rel_perm_w_->Init("wetting", *wrm_list);
  std::string krel_method_name = mp_list_->get<std::string>("relative permeability");
  //rel_perm_w_->ProcessStringRelativePermeability(krel_method_name);

  // Select a proper matrix class. 
  Teuchos::ParameterList& tmp_list = mp_list_->sublist("operators").sublist("diffusion operator");
  Teuchos::ParameterList oplist_matrix = tmp_list.sublist("matrix");
  Teuchos::ParameterList oplist_pc = tmp_list.sublist("preconditioner");

  // create diffusion operators op_matrix_ and op_preconditioner_.
  // They will need to be initialized, which is done later in InitNextTI.
  Operators::OperatorDiffusionFactory opfactory;
  op1_matrix_ = opfactory.Create(oplist_matrix, mesh_, op_bc_p_, rho_, gravity_);
  //op1_matrix_->SetDensity(rho_);
  op1_preconditioner_ = opfactory.Create(oplist_pc, mesh_, op_bc_p_, rho_, gravity_);
  for (int n_op = 0; n_op < soln_->size(); n_op++) {
    ops_.push_back(opfactory.Create(oplist_pc, mesh_, op_bc_p_, rho_, gravity_));
  }
  //op_test_ = opfactory.Create(mesh_, op_bc_p_, oplist_matrix, gravity_, 1);

  Teuchos::ParameterList olist_adv = mp_list_->sublist("operators").sublist("advection operator");
  op2_preconditioner_ = Teuchos::rcp(new Operators::OperatorAdvection(olist_adv, mesh_));

  // accumulation term
  op_acc_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, op2_preconditioner_->global_operator()));

  // create saturation from solution
  saturation_phase1_ = Teuchos::rcp(new CompositeVector(*soln_->SubVector(1)->Data()));
  saturation_phase1_->Scale(-1.0);
  saturation_phase1_->Shift(1.0);

  // init upwind operator
  Teuchos::ParameterList upw_list = mp_list_->sublist("operators")
                                            .sublist("diffusion operator")
                                            .sublist("upwind");
  Operators::UpwindFactory<RelativePermeability> upwind_factory;
  upwind_ = upwind_factory.Create(mesh_, rel_perm_w_, upw_list);

  // init upwind_velocity
  upwind_velocity_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("phase1_velocity", passwd_)));
  tmp_flux_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("phase1_velocity", passwd_)));

} // End InitializePressure()



/* ******************************************************************
* Generic initialization of a next time integration phase.
****************************************************************** */
void Phase1_PK::InitNextTI()
{
  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);
  cvs.AddComponent("face", AmanziMesh::FACE, 1);

  Teuchos::RCP<CompositeVector> Krel = Teuchos::rcp(new CompositeVector(cvs));
  Krel->PutScalarMasterAndGhosted(1.0);

  // Calculate phase relative permeability needed to initialize diffusion operator
  //CompositeVector& phase1_saturation = *S_->GetFieldData("phase1_saturation", passwd_);
  rel_perm_w_->Compute(*saturation_phase1_);
  rel_perm_w_->Krel()->Scale(rho_/mu_);
  rel_perm_w_->dKdS()->Scale(-rho_/mu_);

  upwind_velocity_ = S_->GetFieldData("phase1_velocity", passwd_);
  RelativePermeabilityUpwindFn func1 = &RelativePermeability::Value;
  upwind_->Compute(*upwind_velocity_, *upwind_velocity_, bc_model_s_, bc_value_s_, *rel_perm_w_->Krel(), *rel_perm_w_->Krel(), func1);

  RelativePermeabilityUpwindFn func2 = &RelativePermeability::Derivative;
  upwind_->Compute(*upwind_velocity_, *upwind_velocity_, bc_model_s_, bc_value_s_, *rel_perm_w_->dKdS(), *rel_perm_w_->dKdS(), func2);

  // initialize matrix and preconditioner operators
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_); 
  op1_matrix_->Setup(Kptr, rel_perm_w_->Krel(), Teuchos::null);
  op1_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op1_matrix_->ApplyBCs(true, true); // arguments have no effect for Finite Volume

  op1_preconditioner_->Setup(Kptr, rel_perm_w_->Krel(), Teuchos::null);
  op1_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // Advection
  op2_preconditioner_->Setup(*upwind_velocity_);
  op2_preconditioner_->UpdateMatrices(*rel_perm_w_->dKdS());

  for (int ii = 0; ii < ops_.size(); ii++) {
    ops_[ii]->Setup(Kptr, rel_perm_w_->Krel(), Teuchos::null);
    ops_[ii]->UpdateMatrices(Teuchos::null, Teuchos::null);
  }

} // End InitNextTI

bool Phase1_PK::AdvanceStep(double t_old, double t_new, bool reinit) 
{
  // do nothing
}

/* ******************************************************************* 
* Wrapper for Phase1_PK::CommitState. This is the interface used 
* by the MPC.
******************************************************************* */
void Phase1_PK::CommitStep(double t_old, double t_new)
{
  // Write pressure of phase1 to state
  CompositeVector& p1 = *S_->GetFieldData("phase1_pressure", passwd_);
  p1 = *soln_->SubVector(0)->Data();

  // write saturation of phase1 to state
  CompositeVector& s1 = *S_->GetFieldData("phase1_saturation", passwd_);
  *saturation_phase1_ = *soln_->SubVector(1)->Data();
  saturation_phase1_->Scale(-1.0);
  saturation_phase1_->Shift(1.0);
  s1 = *saturation_phase1_;

  // compute phase1 flux from p1
  CompositeVector& phase1_flux = *S_->GetFieldData("phase1_velocity", passwd_);
  op1_matrix_->UpdateFlux(p1, phase1_flux);
}


/* ******************************************************************
*  Temporary convertion from double to tensor.                                               
****************************************************************** */
void Phase1_PK::SetAbsolutePermeabilityTensor()
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
void Phase1_PK::ComputeBCs() {
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
    bc_value_p_[f] = bc->second;
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
    bc_value_s_[f] = 1.0 - bc->second; // bc value for s1
  }

  for (bc = bc_flux_phase1_->begin(); bc != bc_flux_phase1_->end(); ++bc) {
    int f = bc->first;
    bc_model_p_[f] = Operators::OPERATOR_BC_NEUMANN;
    bc_value_p_[f] = bc->second;
  }

  // mark missing boundary conditions as zero flux conditions only for pressure equation
  AmanziMesh::Entity_ID_List cells;
  missed_bc_faces_ = 0;
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_model_p_[f] == Operators::OPERATOR_BC_NONE) {
      cells.clear();
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
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

}  // namespace Flow
} // End namespace Amanzi
