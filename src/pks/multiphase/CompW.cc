/*
  MultiPhase

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
#include "MFD3D_Diffusion.hh"
#include "OperatorDefs.hh"
#include "CommonDefs.hh"

#include "Mesh.hh"

#include "CompW.hh"
#include "MultiphaseDefs.hh"

// namespace
namespace Amanzi {
namespace Multiphase {

CompW_PK::CompW_PK(Teuchos::ParameterList& pk_tree,
                      const Teuchos::RCP<Teuchos::ParameterList>& specs_list,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& soln):
              bc_pressure_(NULL),
              bc_rhl_(NULL),
              bc_saturation_(NULL),
              src_sink_(NULL),
              S_(S),
              soln_(soln),
              passwd_("state"),
              jacobian_type_("analytic")
{
  comp_list_ = Teuchos::sublist(specs_list, "Component 1");
  if (!specs_list->isSublist("Water retention models")) {
    Errors::Message msg;
    msg << "Component PK: input XML files does not have <Water retention models> sublist in <MPMC Specs> \n";
    Exceptions::amanzi_throw(msg);
  } else {
    wrm_list_ = specs_list->sublist("Water retention models");  
  }

  jacobian_type_ = specs_list->get<std::string>("Jacobian type", "analytic");
  
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

  if (!S->HasField("porosity")) {
    S->RequireScalar("porosity", passwd_);
  }

  if (!S->HasField("gravity")) {
    S->RequireConstantVector("gravity", passwd_, dim_);
  }

  if (!S->HasField("water component density")) {
    S->RequireField("water component density", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField("pressure_w")) {
    S->RequireField("pressure_w", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  
  if (!S->HasField("permeability")) {
    S->RequireField("permeability", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim_);
  }

  if (!S->HasField("diffusion tensor")) {
    S->RequireField("diffustion tensor", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim_);
  }

  if (!S->HasField("saturation_w")) {
    S->RequireField("saturation_w", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField("hydrogen density liquid")) {
    S->RequireField("hydrogen density liquid", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S->HasField("velocity_wet")) {
    S->RequireField("velocity_wet", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      //->SetComponents(names, locations, ndofs);
      ->SetComponent("face", AmanziMesh::FACE, 1);
  }

} // End CompW_PK constructor


/* ******************************************************************
* Clean memory.
****************************************************************** */
CompW_PK::~CompW_PK()
{
  if (bc_pressure_ != NULL) delete bc_pressure_;
  if (bc_saturation_ != NULL) delete bc_saturation_;
  if (bc_rhl_ != NULL) delete bc_rhl_;
  //if (bc_flux_ != NULL) delete bc_flux_;

  if (src_sink_ != NULL) delete src_sink_;
  if (vo_ != NULL) delete vo_;
}


void CompW_PK::InitializeFields()
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

  if (S_->HasField("saturation_w")) {
    if (!S_->GetField("saturation_w", passwd_)->initialized()) {
      S_->GetFieldData("saturation_w", passwd_)->PutScalar(1.0);
      S_->GetField("saturation_w", passwd_)->set_initialized();
    }
  }

  if (S_->HasField("hydrogen density liquid")) {
    if (!S_->GetField("hydrogen density liquid", passwd_)->initialized()) {
      S_->GetFieldData("hydrogen density liquid", passwd_)->PutScalar(0.0);
      S_->GetField("hydrogen density liquid", passwd_)->set_initialized();
    }
  }

  if (!S_->GetField("water component density", passwd_)->initialized()) {
    S_->GetFieldData("water component density", passwd_)->PutScalar(0.0);
    S_->GetField("water component density", passwd_)->set_initialized();
  }

  if (!S_->GetField("pressure_w", passwd_)->initialized()) {
    S_->GetFieldData("pressure_w", passwd_)->PutScalar(0.0);
    S_->GetField("pressure_w", passwd_)->set_initialized();
  }

  if (!S_->GetField("velocity_wet", passwd_)->initialized()) {
    S_->GetFieldData("velocity_wet", passwd_)->PutScalar(0.0);
    S_->GetField("velocity_wet", passwd_)->set_initialized();
  }

}


void CompW_PK::InitializeComponent()
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
  rho_w_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("water component density")));
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
  bc_model_.resize(nfaces_wghost_, 0);
  bc_value_p_.resize(nfaces_wghost_, 0.0);
  bc_value_s_.resize(nfaces_wghost_, 0.0);
  bc_value_rhl_.resize(nfaces_wghost_, 0.0);

  //ProcessParameterList(*comp_list_);
  // Create the BC objects.
  if (!comp_list_->isSublist("boundary conditions")) {
    Errors::Message msg;
    msg << "Comp_PK: problem does not have <boundary conditions> list\n";
    Exceptions::amanzi_throw(msg);
  }
  Teuchos::RCP<Teuchos::ParameterList>
        bc_list = Teuchos::rcp(new Teuchos::ParameterList(comp_list_->sublist("boundary conditions", true)));
  /* FIXME
  MultiphaseBCFactory bc_factory(mesh_, bc_list);
  bc_rhl_ = bc_factory.CreateHydrogenDensity(bc_submodel_);
  bc_pressure_ = bc_factory.CreatePressure(bc_submodel_);
  bc_saturation_ = bc_factory.CreateSaturation(bc_submodel_);
  */

  //double time = S_->time();
  //if (time >= 0.0) T_physics_ = time;
  double time = 0.0;

  time = T_physics_;
  bc_pressure_->Compute(time, time);
  bc_saturation_->Compute(time, time);
  bc_rhl_->Compute(time, time);

  // Create water retention models.
  rel_perm_w_ = Teuchos::rcp(new MPCoeff(mesh_));
  rel_perm_w_->Init("wetting", wrm_list_);

  // rel_perm_w_->ProcessStringRelativePermeability(krel_method_name);

  ComputeBCs();

  // create operator boundary condition object
  op_bc_p_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  op_bc_s_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  op_bc_rhl_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

  // allocate memory for absolute permeability
  K_.resize(ncells_wghost_);
  D1_.resize(ncells_wghost_);
  SetAbsolutePermeabilityTensor();
  SetDiffusionTensor();

  // if (src_sink_distribution_ & Amanzi::CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
  //   Kxy_ = Teuchos::rcp(new Epetra_Vector(mesh_->cell_map(false)));
  // }

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
  op2_matrix_ = opfactory.Create(oplist_matrix, mesh_, op_bc_rhl_, 1.0, gravity_);

  if (jacobian_type_ == "numerical")
    for (int n_op = 0; n_op < soln_->size(); n_op++) {
      ops_.push_back(opfactory.Create(oplist_pc, mesh_, op_bc_p_, 1.0, gravity_));
    }

  op1_preconditioner_ = opfactory.Create(oplist_pc, mesh_, op_bc_p_, 1.0, gravity_);
  op3_preconditioner_ = opfactory.Create(oplist_pc, mesh_, op_bc_rhl_, 1.0, gravity_);
  Teuchos::ParameterList olist_adv = comp_list_->sublist("operators").sublist("advection operator");
  op2_preconditioner_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, mesh_));
  op_prec_sat_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(olist_adv, op2_preconditioner_->global_operator()));

  op_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op_prec_sat_->global_operator()));

  // initialize saturation of phase1 to compute flux of phase1 later
  //saturation_phase1_ = Teuchos::rcp(new CompositeVector(soln_->SubVector(1)->Data()->Map()));

  // init upwind operator
  Teuchos::ParameterList upw_list = comp_list_->sublist("operators")
                                            .sublist("diffusion operator")
                                            .sublist("upwind");

  //Operators::UpwindFactory<RelativePermeability> upwind_factory;
  //upwind_n_ = upwind_factory.Create(mesh_, rel_perm_w_, upw_list);

  upwind_w_ = Teuchos::rcp(new Operators::UpwindFlux<MPCoeff>(mesh_, rel_perm_w_));
  upwind_w_->Init(upw_list);
  //std::cout << "Done creating upwind operator \n";

  // init upwind_velocity
  upwind_vw_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("velocity_wet", passwd_)));
  //std::cout << "Done init upwind_velocity \n";

  tmp_flux_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("velocity_wet", passwd_)));
  //std::cout << "Done InitializeComponent\n";
} 


/* ******************************************************************
* Generic initialization of a next time integration phase.
****************************************************************** */
void CompW_PK::InitNextTI()
{
  // initialize operators
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > D1ptr = Teuchos::rcpFromRef(D1_);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->Setup(Kptr, rel_perm_w_->Krel(), Teuchos::null);
  op1_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);

  op2_matrix_->Setup(D1ptr, Teuchos::null, Teuchos::null);
  op2_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);

  if (jacobian_type_ == "numerical")
    for (int ii = 0; ii < ops_.size(); ii++) {
      ops_[ii]->Setup(Kptr, rel_perm_w_->Krel(), Teuchos::null);
      ops_[ii]->UpdateMatrices(Teuchos::null, Teuchos::null);
    }
  //op1_preconditioner_->Setup(Kptr, rel_perm_w_->Krel(), Teuchos::null);
  //op1_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // Advection
  //op2_preconditioner_->Setup(*upwind_vw_);
  //op2_preconditioner_->UpdateMatrices(*rel_perm_w_->dKdS());

  // diffusion hydrogen
  //op3_preconditioner_->Setup(Kptr, Teuchos::null, Teuchos::null);
  //op3_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);
}


/* ******************************************************************* 
* This is the interface used by the MPC.
******************************************************************* */
void CompW_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S)
{
  // Write pressure to state
  auto p1 = soln_->SubVector(0)->Data();

  // write saturation of phase1 to state
  CompositeVector& s1 = *S_->GetFieldData("saturation_w", passwd_);
  s1 = *soln_->SubVector(1)->Data();
  //s1.Print(std::cout);

  // compute upwind velocities from pressure 
  // Calculate total mobility needed to initialize diffusion operator
  rel_perm_w_->Compute(s1);
  //std::cout << "krel_w before upwind " << *rel_perm_w_->Krel()->ViewComponent("cell") << "\n";
  upwind_vw_ = S_->GetFieldData("velocity_wet", passwd_);
  //std::cout << "upwind velocity: " << *upwind_vw_->ViewComponent("face") << "\n";

  upwind_w_->Compute(*upwind_vw_, *upwind_vw_, bc_model_, *rel_perm_w_->Krel());
  rel_perm_w_->Krel()->Scale((t_new - t_old)/mu_); // don't need rho_ when compute darcy flux

  auto phase1_flux = S_->GetFieldData("velocity_wet", passwd_);
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->global_operator()->Init();
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->Setup(Kptr, rel_perm_w_->Krel(), Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->UpdateMatrices(Teuchos::null, Teuchos::null);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->ApplyBCs(true, true ,true);
  ((Teuchos::RCP<Operators::PDE_Diffusion>)op1_matrix_)->UpdateFlux(p1.ptr(), phase1_flux.ptr());
}


/* ******************************************************************
*  Temporary convertion from double to tensor.                                               
****************************************************************** */
void CompW_PK::SetAbsolutePermeabilityTensor()
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
void CompW_PK::SetDiffusionTensor()
{
  const CompositeVector& cv1 = *S_->GetFieldData("diffusion tensor");
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
}


// need to fix flag essential bc. We need Dirichlet bc for all 4 primary
// variables before incrementing Dirichlet bc counter 
/* ******************************************************************* 
* Compute boundary conditions based on bc type and values
******************************************************************* */
void CompW_PK::ComputeBCs() {
  int flag_essential_bc = 0;
  dirichlet_bc_faces_ = 0;

  for (int n = 0; n < bc_model_.size(); n++) {
    bc_model_[n] = Operators::OPERATOR_BC_NONE;
    bc_value_p_[n] = 0.0;
    bc_value_s_[n] = 0.0;
    bc_value_rhl_[n] = 0.0;
  }

  Flow::FlowBoundaryFunction::Iterator bc;
  for (bc = bc_pressure_->begin(); bc != bc_pressure_->end(); ++bc) {
    int f = bc->first;
    bc_model_[f] = Operators::OPERATOR_BC_DIRICHLET;
    bc_value_p_[f] = bc->second[0];
    flag_essential_bc = 1;
    dirichlet_bc_faces_++;
  }

  for (bc = bc_saturation_->begin(); bc != bc_saturation_->end(); ++bc) {
    int f = bc->first;
    bc_value_s_[f] = bc->second[0];
  }

  for (bc = bc_rhl_->begin(); bc != bc_rhl_->end(); ++bc) {
    int f = bc->first;
    bc_value_rhl_[f] = bc->second[0];
  }

  // mark missing boundary conditions as zero flux conditions 
  AmanziMesh::Entity_ID_List cells;
  missed_bc_faces_ = 0;
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_model_[f] == Operators::OPERATOR_BC_NONE) {
      cells.clear();
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int ncells = cells.size();

      if (ncells == 1) {
        bc_model_[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value_p_[f] = 0.0;
        bc_value_s_[f] = 0.0;
        bc_value_rhl_[f] = 0.0;
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


void CompW_PK::DeriveFaceValuesFromCellValues(
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

}  // namespace Flow
} // End namespace Amanzi
