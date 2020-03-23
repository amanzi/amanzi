/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Daniil Svyatsiyy
          Konstantin Lipnikov (lipnikov@lanl.gov)

  Major transport algorithms.
*/

#include <algorithm>
#include <vector>

#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "BCs.hh"
#include "errors.hh"
#include "FieldEvaluator.hh"
#include "LinearOperatorDefs.hh"
#include "LinearOperatorFactory.hh"
#include "Mesh.hh"
#include "PDE_Accumulation.hh"
#include "PDE_AdvectionUpwindFactory.hh"
#include "PDE_Diffusion.hh"
#include "PDE_DiffusionFactory.hh"
#include "PK_DomainFunctionFactory.hh"
#include "PK_Utils.hh"
#include "WhetStoneDefs.hh"

// Amanzi::Transport
#include "MultiscaleTransportPorosityFactory.hh"
#include "TransportImplicit_PK.hh"
#include "TransportBoundaryFunction_Alquimia.hh"
#include "TransportDomainFunction.hh"
#include "TransportSourceFunction_Alquimia.hh"

namespace Amanzi {
namespace Transport {

/* ******************************************************************
* New constructor compatible with new MPC framework.
****************************************************************** */
TransportImplicit_PK::TransportImplicit_PK(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& glist,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& soln) :
    Transport_PK(pk_tree, glist, S, soln)
{
  if (tp_list_->isSublist("time integrator"))
    ti_list_ = Teuchos::sublist(tp_list_, "time integrator", true);

  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "solvers", true);
}


/* ******************************************************************
* Simple constructor for unit tests.
****************************************************************** */
TransportImplicit_PK::TransportImplicit_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                           Teuchos::RCP<State> S, 
                                           const std::string& pk_list_name,
                                           std::vector<std::string>& component_names) :
  Transport_PK(glist, S, pk_list_name, component_names)
{

  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  tp_list_ = Teuchos::sublist(pk_list, pk_list_name, true);

  // We also need miscaleneous sublists
  preconditioner_list_ = Teuchos::sublist(glist, "preconditioners", true);
  linear_operator_list_ = Teuchos::sublist(glist, "solvers", true);
  if (tp_list_->isSublist("time integrator"))
    ti_list_ = Teuchos::sublist(tp_list_, "time integrator", true);
}

  
/* ******************************************************************
* Initialization
****************************************************************** */
void TransportImplicit_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  Transport_PK::Initialize(S);  

  // domain name
  Key domain = tp_list_->template get<std::string>("domain name", "domain");
  auto vo_list = tp_list_->sublist("verbose object"); 
  vo_ = Teuchos::rcp(new VerboseObject("TransportImpl-" + domain, *tp_list_)); 

  // Create pointers to the primary solution field pressure.
  const auto& solution = S->GetFieldData(tcc_key_, "state");
  soln_->SetData(solution); 
  
  // boundary conditions
  op_bc_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  auto& values = op_bc_->bc_value();
  auto& models = op_bc_->bc_model();

  Teuchos::ParameterList& oplist = tp_list_->sublist("operators")
                                            .sublist("advection operator")
                                            .sublist("matrix");

  Operators::PDE_AdvectionUpwindFactory adv_factory;
  op_adv_ = adv_factory.Create(oplist, mesh_);

  op_ = op_adv_->global_operator();
  op_adv_->SetBCs(op_bc_, op_bc_);

  // refresh data BC and source data  
  UpdateBoundaryData(t_physics_, t_physics_, 0);

  auto flux = S->GetFieldData(darcy_flux_key_);
  op_adv_->Setup(*flux);
  op_adv_->UpdateMatrices(flux.ptr());

  Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh_)->SetGhosted(true)->AddComponent("cell", AmanziMesh::CELL, 1);

  acc_term_ = Teuchos::rcp(new CompositeVector(*cvs));
  acc_term_prev_ = Teuchos::rcp(new CompositeVector(*cvs));
  op_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op_));

  op_->SymbolicAssembleMatrix();
  op_->CreateCheckPoint();
  // op_->AssembleMatrix(); 

  // generic linear solver
  if (ti_list_ != Teuchos::null) {
    // preconditioner
    if (ti_list_->isParameter("linear solver"))
      solver_name_ = ti_list_->get<std::string>("linear solver");

    if (ti_list_->isParameter("preconditioner")) {
      std::string name = ti_list_->get<std::string>("preconditioner");
      Teuchos::ParameterList pc_list = preconditioner_list_->sublist(name);
      op_->InitializePreconditioner(pc_list);
    }
  }
  
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << vo_->color("green") << "Initialization of PK is complete." 
               << vo_->reset() << std::endl << std::endl;
  }
}


/* ******************************************************************* 
* Performs one time step from t_old to t_new. The boundary conditions
* are calculated only once, during the initialization step.  
******************************************************************* */
bool TransportImplicit_PK::AdvanceStep(double t_old, double t_new, bool reinit) 
{
  bool fail = false;
  dt_ = t_new - t_old;
  double dt_MPC(dt_);
  
  // populating next state of concentrations
  tcc->ScatterMasterToGhosted("cell");

  CompositeVector tcc_aux(*acc_term_);

  int icount(0);
  std::string name;

  for (int i = 0; i < num_aqueous; i++) {
    op_->RestoreCheckPoint();
    // op_acc_->local_op(0)->Rescale(0.0);
  
    // add accumulation term
    Epetra_MultiVector& acc_term_c = *acc_term_->ViewComponent("cell");
    Epetra_MultiVector& acc_term_prev_c = *acc_term_prev_->ViewComponent("cell");
  
    for (int c = 0; c < ncells_owned; c++) {
      acc_term_c[0][c] = (*phi)[0][c] * (*ws)[0][c];
      acc_term_prev_c[0][c] = (*phi)[0][c] * (*ws_prev)[0][c];    
    }

    *(*tcc_aux.ViewComponent("cell"))(0) = *(*tcc->ViewComponent("cell"))(i);  
    op_acc_->AddAccumulationDelta(tcc_aux, *acc_term_prev_, *acc_term_, dt_, "cell");

    // refresh data BC and source data  
    UpdateBoundaryData(t_old, t_new, i);

    CompositeVector& rhs = *op_->rhs();
    Epetra_MultiVector& rhs_cell = *rhs.ViewComponent("cell");
  
    // apply boundary conditions
    op_adv_->UpdateMatrices(S_->GetFieldData(darcy_flux_key_).ptr());
    op_adv_->ApplyBCs(true, true, true);

    // add sources
    ComputeSources_(t_new, t_new - t_old, rhs_cell, *tcc->ViewComponent("cell"), i, i);
  
    // assemble matrix
    op_->AssembleMatrix();
    op_->UpdatePreconditioner();

    // create linear solver and calculate new pressure
    AmanziSolvers::LinearOperatorFactory<Operators::Operator, CompositeVector, CompositeVectorSpace> factory;
    Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::Operator, CompositeVector, CompositeVectorSpace> >
       solver = factory.Create(solver_name_, *linear_operator_list_, op_);

    solver->add_criteria(AmanziSolvers::LIN_SOLVER_MAKE_ONE_ITERATION);
    solver->ApplyInverse(rhs, tcc_aux);

    *(*tcc_tmp->ViewComponent("cell"))(i) = *(*tcc_aux.ViewComponent("cell"))(0);

    // statistics
    icount += solver->num_itrs();
    name = solver->name();
  }

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    double sol_norm;
    tcc_tmp->Norm2(&sol_norm);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "transport solver (" << name
               << "): ||sol||=" << sol_norm 
               << "  avg itrs=" << icount / num_aqueous << std::endl;
      VV_PrintSoluteExtrema(*tcc->ViewComponent("cell"), t_new - t_old, "");
  }

  // estimate time multiplier
  // dt_desirable_ = ts_control_->get_timestep(dt_MPC, 1);

  return fail;
}


/* ******************************************************************* 
* Boundary data
******************************************************************* */
void TransportImplicit_PK::UpdateBoundaryData(double t_old, double t_new, int component)
{
  for (int i = 0; i < bcs_.size(); i++) {
    bcs_[i]->Compute(t_old, t_new);
  }

  auto& values = op_bc_->bc_value();
  auto& models = op_bc_->bc_model();

  ComputeBCs_(models, values, component);
}

}  // namespace Transport
}  // namespace Amazni
