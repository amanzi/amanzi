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
#include "Mesh.hh"

#include "Pressure_PK.hh"

// namespace
namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Simple construction 
****************************************************************** */
Pressure_PK::Pressure_PK(Teuchos::ParameterList& pk_tree,
                      const Teuchos::RCP<Teuchos::ParameterList>& glist,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& soln):
              bc_pressure(NULL),
              bc_flux(NULL),
              vo_(NULL),
              S_(S),
              sol_tree_(soln),
              rho_(1.0), mu_(1.0),
              passwd_("state")
{
  Teuchos::RCP<Teuchos::ParameterList> flow_list = Teuchos::sublist(glist, "Flow", true);
  mp_list_ = Teuchos::sublist(flow_list, "Pressure Equation", true);
  linear_operator_list_ = Teuchos::sublist(glist, "Solvers", true);
  preconditioner_list_ = Teuchos::sublist(glist, "Preconditioners", true);
}

void Pressure_PK::Setup()
{
  mesh_ = S_->GetMesh();
  dim_ = mesh_->space_dimension();
  std::vector<std::string> names(2);
  names[0] = "cell";
  names[1] = "face";

  std::vector<int> ndofs(2, 1);

  std::vector<AmanziMesh::Entity_kind> locations(2);
  locations[0] = AmanziMesh::CELL;
  locations[1] = AmanziMesh::FACE;

  // require state variables for the Multiphase Pressure
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
  
  if (!S_->HasField("gravity")) {
    S_->RequireConstantVector("gravity", passwd_, dim_);
  }

  if (!S_->HasField("pressure")) {
    S_->RequireField("pressure", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  
  if (!S_->HasField("permeability")) {
    S_->RequireField("permeability", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, dim_);
  }

  if (!S_->HasField("water_saturation")) {
    S_->RequireField("water_saturation", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  }

  // darcy_flux is defined on faces, but we also define it on cells to make it
  // compatible with fractional flow for multiplication.
  if (!S_->HasField("darcy_flux")) {
    S_->RequireField("darcy_flux", passwd_)->SetMesh(mesh_)->SetGhosted(true)
      ->SetComponents(names, locations, ndofs);
  }
}


/* ******************************************************************
* Clean memory.
****************************************************************** */
Pressure_PK::~Pressure_PK()
{
  if (bc_pressure != NULL) delete bc_pressure;
  if (bc_flux != NULL) delete bc_flux;

  for (int i = 0; i < srcs_.size(); i++) {
    if (srcs_[i] != NULL) delete srcs_[i]; 
  }
  if (vo_ != NULL) delete vo_;
}


/* ******************************************************************
* Initialization of fields
****************************************************************** */
void Pressure_PK::InitializeFields()
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

  if (S_->HasField("water_saturation")) {
    if (!S_->GetField("water_saturation", passwd_)->initialized()) {
      S_->GetFieldData("water_saturation", passwd_)->PutScalar(0.0);
      S_->GetField("water_saturation", passwd_)->set_initialized();
    }
  }

  if (!S_->GetField("pressure", passwd_)->initialized()) {
    S_->GetFieldData("pressure", passwd_)->PutScalar(0.0);
    S_->GetField("pressure", passwd_)->set_initialized();
  }

  if (!S_->GetField("darcy_flux", passwd_)->initialized()) {
    S_->GetFieldData("darcy_flux", passwd_)->PutScalar(0.0);
    S_->GetField("darcy_flux", passwd_)->set_initialized();
  }
}


void Pressure_PK::InitializePressure()
{
  // times
  dt_ = -1.0;

  // Initilize various common data depending on mesh and state.
  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);

  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);

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

  // create verbosity object
  vo_ = new VerboseObject("MultiphasePK::Pressure", *mp_list_); 

  // Allocate memory for boundary data.
  bc_submodel.resize(nfaces_wghost, 0);
  
  ProcessParameterList(*mp_list_);

  // Create the source object if any
  /*
  if (mp_list_->isSublist("source terms")) {
    Teuchos::RCP<Teuchos::ParameterList> src_list = Teuchos::rcpFromRef(mp_list_.sublist("source terms", true));
    FlowSourceFactory factory(mesh_, src_list);
    factory.Create(srcs_);
  }
  */

  // compute boundary values, which can be functions of time
  // but we only use boundary conditions independent of time for now
  bc_pressure->Compute(0.0, 0.0);
  bc_flux->Compute(0.0, 0.0);
  ComputeBCs();

  // create operator boundary condition object
  op_bc_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

  // allocate memory for absolute permeability
  K_.resize(ncells_wghost);

  // Create water retention models.
  if (!mp_list_->isSublist("Water retention models")) {
    Errors::Message msg;
    msg << "Multiphase PK: input XML files does not have <Water retention models> sublist in <Multiphase problem> \n";
    Exceptions::amanzi_throw(msg);
  }
  Teuchos::RCP<Teuchos::ParameterList> wrm_list = Teuchos::sublist(mp_list_, "Water retention models", true);
  rel_perm_w_ = Teuchos::rcp(new RelativePermeability(mesh_));
  rel_perm_n_ = Teuchos::rcp(new RelativePermeability(mesh_));
  rel_perm_w_->Init("wetting", *wrm_list);
  rel_perm_n_->Init("non wetting", *wrm_list);

  std::string krel_method_name = mp_list_->get<std::string>("relative permeability");
  //rel_perm_w_->ProcessStringRelativePermeability(krel_method_name);
  //rel_perm_n_->ProcessStringRelativePermeability(krel_method_name);

  // initialize total mobility from relative permeability
  tot_mobility_ = Teuchos::rcp(new CompositeVector(rel_perm_w_->Krel()->Map()));
  tot_mobility_->PutScalar(0.0);

  // Select a proper matrix class. 
  Teuchos::ParameterList& tmp_list = mp_list_->sublist("operators").sublist("diffusion operator");
  Teuchos::ParameterList& oplist_matrix = tmp_list.sublist("matrix");
  Teuchos::ParameterList& oplist_pc = tmp_list.sublist("preconditioner");

  // create diffusion operators op_matrix_ and op_preconditioner_.
  // They will need to be initialized, which is done later in InitNextTI.
  Operators::PDE_DiffusionFactory opfactory;
  op_matrix_ = opfactory.Create(oplist_matrix, mesh_, op_bc_, 1.0, gravity_);
  op_preconditioner_ = Teuchos::rcp(new Operators::PDE_DiffusionFV(oplist_pc, mesh_));

  // Create the solution (pressure) vector and auxiliary vector for time history.
  sol_ = Teuchos::rcp(new CompositeVector(*S_->GetFieldData("pressure",passwd_)));
  sol_tree_->SetData(sol_);

  // Calculate total mobility needed to initialize diffusion operator
  CompositeVector water_saturation = *S_->GetFieldData("water_saturation");
  rel_perm_w_->Compute(water_saturation);
  rel_perm_n_->Compute(water_saturation);
  *tot_mobility_ = *rel_perm_w_->Krel();
  tot_mobility_->Update(1.0/mu_[1], *rel_perm_n_->Krel(), 1.0/mu_[0]);

  gravity_factor_ = Teuchos::rcp(new CompositeVector(*tot_mobility_));
  include_gravity_ = oplist_pc.get<bool>("gravity", false);
  //if (include_gravity_) op_preconditioner_->SetGravity(gravity_);
  op_preconditioner_->SetBCs(op_bc_, op_bc_);
  
  // We want to compute the total mobility lambda_tot on the faces using 
  // harmonic average. Since Krel is upwinded by OperatorDiffusion by default, 
  // and absolute permeability K is computed using harmonic average, we incorporate 
  // total mobility into absolute permeability.
  SetAbsolutePermeabilityTensor(tot_mobility_);

  // initialize matrix and preconditioner operators
  Teuchos::RCP<std::vector<WhetStone::Tensor> > Kptr = Teuchos::rcpFromRef(K_); 
  op_matrix_->Setup(Kptr, Teuchos::null, Teuchos::null);
  op_matrix_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op_matrix_->ApplyBCs(true, true, true);

  op_preconditioner_->Setup(Kptr, Teuchos::null, Teuchos::null);
  op_preconditioner_->UpdateMatrices(Teuchos::null, Teuchos::null);
  op_preconditioner_->ApplyBCs(true, true, true);
  op_ = op_preconditioner_->global_operator();

  // generic linear solver
  solver_name_ = mp_list_->get<std::string>("linear solver");

  // preconditioner
  preconditioner_name_ = mp_list_->get<std::string>("preconditioner");
  ASSERT(preconditioner_list_->isSublist(preconditioner_name_));

  // Initialize source
  /*
  if (src_sink_ != NULL) {
    if (src_sink_distribution_ & Amanzi::CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      src_sink_->ComputeDistribute(t_old, t_new, Kxy->Values()); 
    } else {
      src_sink_->ComputeDistribute(t_old, t_new);
    }
  }
  */
} // End InitializePressure()


/* ******************************************************************* 
* Solve the pressure equation
******************************************************************* */
bool Pressure_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  double dt_ = t_new - t_old;

  // Update source and boundary conditions, should have no effect in this case, since pressure
  // is steady state
  /*
  if (src_sink_ != NULL) {
    if (src_sink_distribution_ & Amanzi::CommonDefs::DOMAIN_FUNCTION_ACTION_DISTRIBUTE_PERMEABILITY) {
      src_sink_->ComputeDistribute(t_old, t_new, Kxy->Values()); 
    } else {
      src_sink_->ComputeDistribute(t_old, t_new);
    }
  }
  */

  bc_pressure->Compute(t_old, t_new);
  bc_flux->Compute(t_old, t_new);
  ComputeBCs();

  // Update the preconditioner and solve the problem
  UpdatePreconditioner(1.0, sol_tree_, 1.0);
  Teuchos::RCP<TreeVector> rhs = Teuchos::rcp(new TreeVector());
  // rhs_ is computed from op_preconditioner_
  rhs->SetData(rhs_);
  ApplyPreconditioner(rhs, sol_tree_);
  
  return false;
}


/* ******************************************************************* 
* Write back to solution and other data to the state
******************************************************************* */
void Pressure_PK::CommitState(const Teuchos::Ptr<State>& S)
{
  // update solution
  auto p = S->GetFieldData("pressure", passwd_);
  *p = *sol_tree_->Data();

  // calculate and update darcy mass flux
  auto darcy_flux = S->GetFieldData("darcy_flux", passwd_);
  op_preconditioner_->UpdateFlux(p.ptr(), darcy_flux.ptr());

  // usually the flux on the faces needs to be divided by density rho, 
  // but not in pressure equation since we divide both sides by rho
  // in the decoupled formulation
}


/* ******************************************************************* 
* Wrapper for Pressure_PK::CommitState. This is the interface used 
* by the MPC.
******************************************************************* */
void Pressure_PK::CommitStep(double t_old, double t_new)
{
  this->CommitState(S_.ptr());
}


/* ******************************************************************
* Add source and sink terms.                                   
****************************************************************** */
void Pressure_PK::AddSourceTerms(CompositeVector& rhs)
{
  // Do nothing for now
  /*
  Epetra_MultiVector& rhs_cell = *rhs.ViewComponent("cell");
  Flow::FlowDomainFunction::Iterator src;
  for (src = src_sink_->begin(); src != src_sink_->end(); ++src) {
    int c = src->first;
    rhs_cell[0][c] += mesh_->cell_volume(c) * src->second;
  }
  */
}


/* ******************************************************************
*  Temporary convertion from double to tensor.                                               
****************************************************************** */
void Pressure_PK::SetAbsolutePermeabilityTensor(Teuchos::RCP<CompositeVector> lambda)
{
  const CompositeVector& cv = *S_->GetFieldData("permeability");
  cv.ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& perm = *cv.ViewComponent("cell", true);
  const Epetra_MultiVector& lambda_cell = *(lambda->ViewComponent("cell", true));

  if (dim_ == 2) {
    for (int c = 0; c < K_.size(); c++) {
      if (perm[0][c] == perm[1][c]) {
        K_[c].Init(dim_, 1);
        K_[c](0, 0) = perm[0][c] * lambda_cell[0][c];
      } else {
        K_[c].Init(dim_, 2);
        K_[c](0, 0) = perm[0][c] * lambda_cell[0][c];
        K_[c](1, 1) = perm[1][c] * lambda_cell[0][c];
      }
    }    
  } else if (dim_ == 3) {
    for (int c = 0; c < K_.size(); c++) {
      if (perm[0][c] == perm[1][c] && perm[0][c] == perm[2][c]) {
        K_[c].Init(dim_, 1);
        K_[c](0, 0) = perm[0][c] * lambda_cell[0][c];
      } else {
        K_[c].Init(dim_, 2);
        K_[c](0, 0) = perm[0][c] * lambda_cell[0][c];
        K_[c](1, 1) = perm[1][c] * lambda_cell[0][c];
        K_[c](2, 2) = perm[2][c] * lambda_cell[0][c];
      }
    }        
  }
}


/* ******************************************************************* 
* Compute boundary conditions based on bc type and values
******************************************************************* */
void Pressure_PK::ComputeBCs()
{
  std::vector<int>& bc_model = op_bc_->bc_model();
  std::vector<double>& bc_value = op_bc_->bc_value();
  std::vector<double>& bc_mixed = op_bc_->bc_mixed();

  int flag_essential_bc = 0;
  dirichlet_bc_faces_ = 0;

  for (int n = 0; n < bc_model.size(); n++) {
    bc_model[n] = Operators::OPERATOR_BC_NONE;
    bc_value[n] = 0.0;
  }

  Flow::FlowBoundaryFunction::Iterator bc;
  for (bc = bc_pressure->begin(); bc != bc_pressure->end(); ++bc) {
    int f = bc->first;
    bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
    bc_value[f] = bc->second[0];
    flag_essential_bc = 1;
    dirichlet_bc_faces_++;
  }

  for (bc = bc_flux->begin(); bc != bc_flux->end(); ++bc) {
    int f = bc->first;
    bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
    bc_value[f] = bc->second[0];
  }

  // mark missing boundary conditions as zero flux conditions only for pressure equation
  AmanziMesh::Entity_ID_List cells;
  missed_bc_faces_ = 0;
  int nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_model[f] == Operators::OPERATOR_BC_NONE) {
      cells.clear();
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int ncells = cells.size();

      if (ncells == 1) {
        bc_model[f] = Operators::OPERATOR_BC_NEUMANN;
        bc_value[f] = 0.0;
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
