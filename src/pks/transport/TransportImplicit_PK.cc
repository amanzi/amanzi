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
#include "Mesh.hh"
#include "PDE_Accumulation.hh"
#include "PDE_AdvectionUpwindFactory.hh"
#include "PDE_Diffusion.hh"
#include "PDE_DiffusionFactory.hh"
#include "PK_DomainFunctionFactory.hh"
#include "PK_Utils.hh"
#include "WhetStoneDefs.hh"
#include "InverseFactory.hh"

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
void TransportImplicit_PK::Initialize()
{
  Transport_PK::Initialize();  

  // domain name
  Key domain = tp_list_->template get<std::string>("domain name", "domain");
  auto vo_list = tp_list_->sublist("verbose object"); 
  vo_ = Teuchos::rcp(new VerboseObject("TransportImpl-" + domain, *tp_list_)); 

  // Create pointers to the primary solution field tcc.
  auto solution = S_->GetPtrW<CompositeVector>(tcc_key_, Tags::DEFAULT, "state");
  soln_->SetData(solution); 
  
  // boundary conditions
  op_bc_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));
  op_bc_->bc_value();  // allocate internal
  op_bc_->bc_model();  // memory

  // operators
  // -- dispertion and/or diffusion
  if (use_dispersion_) {
    D_.resize(ncells_owned);
    Teuchos::RCP<std::vector<WhetStone::Tensor> > Dptr = Teuchos::rcpFromRef(D_);
    Teuchos::ParameterList& oplist_d = tp_list_->sublist("operators")
                                                .sublist("diffusion operator")
                                                .sublist("matrix");

    Operators::PDE_DiffusionFactory diff_factory;
    op_diff_ = diff_factory.Create(oplist_d, mesh_);
    op_diff_->Setup(Dptr, Teuchos::null, Teuchos::null);
    op_diff_->SetBCs(op_bc_, op_bc_);
    op_ = op_diff_->global_operator();
  }

  // -- advection
  Teuchos::ParameterList& oplist_a = tp_list_->sublist("operators")
                                              .sublist("advection operator")
                                              .sublist("matrix");

  Operators::PDE_AdvectionUpwindFactory adv_factory;
  if (op_.get()) {
    op_adv_ = adv_factory.Create(oplist_a, op_);
  } else {
    op_adv_ = adv_factory.Create(oplist_a, mesh_);
    op_ = op_adv_->global_operator();
  }
  op_adv_->SetBCs(op_bc_, op_bc_);

  // refresh data BC and source data  
  UpdateBoundaryData(t_physics_, t_physics_, 0);

  auto flux = S_->GetPtr<CompositeVector>(vol_flowrate_key_, Tags::DEFAULT);
  op_adv_->Setup(*flux);
  op_adv_->UpdateMatrices(flux.ptr());

  op_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op_));

  op_->CreateCheckPoint();

  // generic linear solver or preconditioner
  if (ti_list_ != Teuchos::null) {
    solver_name_ = ti_list_->get<std::string>("linear solver", "");
    auto pc_name = ti_list_->get<std::string>("preconditioner", "");
    auto inv_list = AmanziSolvers::mergePreconditionerSolverLists(
        pc_name, *preconditioner_list_,
        solver_name_, *linear_solver_list_, true);
    op_->set_inverse_parameters(inv_list);
    op_->InitializeInverse();
  }
  
  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "matrix:" << std::endl
               << op_->PrintDiagnostics() << std::endl
               << vo_->color("green") << "Initialization of PK is complete." 
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
  
  // populating next state of concentrations
  tcc->ScatterMasterToGhosted("cell");

  S_->GetEvaluator(water_content_key_).Update(*S_, "transport");
  const auto& wc = S_->Get<CompositeVector>(water_content_key_, Tags::DEFAULT);
  const auto& wc_prev = S_->Get<CompositeVector>(prev_water_content_key_, Tags::DEFAULT);

  const auto& wc_c = *wc.ViewComponent("cell");
  const auto& sat_c = *S_->Get<CompositeVector>(saturation_liquid_key_, Tags::DEFAULT).ViewComponent("cell");

  int icount(0);
  CompositeVector tcc_aux(wc);

  for (int i = 0; i < num_aqueous; i++) {
    op_->RestoreCheckPoint();
    // op_acc_->local_op(0)->Rescale(0.0);
  
    *(*tcc_aux.ViewComponent("cell"))(0) = *(*tcc->ViewComponent("cell"))(i);  
    op_acc_->AddAccumulationDelta(tcc_aux, wc_prev, wc, dt_, "cell");

    // refresh data BC and source data  
    UpdateBoundaryData(t_old, t_new, i);

    CompositeVector& rhs = *op_->rhs();
    Epetra_MultiVector& rhs_cell = *rhs.ViewComponent("cell");
  
    // apply boundary conditions
    op_adv_->UpdateMatrices(S_->GetPtr<CompositeVector>(vol_flowrate_key_, Tags::DEFAULT).ptr());
    op_adv_->ApplyBCs(true, true, true);

    if (use_dispersion_) {
      int phase;
      double md;
      CalculateDispersionTensor_(*transport_phi, wc_c);
      FindDiffusionValue(component_names_[i], &md, &phase);
      if (md != 0.0) CalculateDiffusionTensor_(md, phase, *transport_phi, sat_c);

      op_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
      op_diff_->ApplyBCs(true, true, true);
    }

    // add sources
    ComputeSources_(t_new, t_new - t_old, rhs_cell, *tcc->ViewComponent("cell"), i, i);
  
    op_->ComputeInverse();
    op_->ApplyInverse(rhs, tcc_aux);

    *(*tcc_tmp->ViewComponent("cell"))(i) = *(*tcc_aux.ViewComponent("cell"))(0);

    // statistics
    icount += op_->num_itrs();
  }

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    double sol_norm;
    tcc_tmp->Norm2(&sol_norm);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "transport solver (" << solver_name_
               << "): ||sol||=" << sol_norm 
               << "  avg itrs=" << icount / num_aqueous << std::endl;
      VV_PrintSoluteExtrema(*tcc->ViewComponent("cell"), t_new - t_old, "");
  }

  // estimate time multiplier
  // dt_desirable_ = ts_control_->get_timestep(dt_MPC, 1);

  return fail;
}


/* ******************************************************************* 
* Recompute linear system for the given component.
******************************************************************* */
void TransportImplicit_PK::UpdateLinearSystem(double t_old, double t_new, int component) 
{
  S_->GetEvaluator(water_content_key_).Update(*S_, "transport");
  const auto& wc = S_->Get<CompositeVector>(water_content_key_, Tags::DEFAULT);
  const auto& wc_prev = S_->Get<CompositeVector>(prev_water_content_key_, Tags::DEFAULT);

  const auto& wc_c = *wc.ViewComponent("cell");
  const auto& sat_c = *S_->Get<CompositeVector>(saturation_liquid_key_, Tags::DEFAULT).ViewComponent("cell");

  CompositeVector tcc_aux(wc);
  *(*tcc_aux.ViewComponent("cell"))(0) = *(*tcc->ViewComponent("cell"))(component);

  double dt = t_new - t_old;
  op_->RestoreCheckPoint();

  op_->rhs()->PutScalar(0.0);
  // op_acc_->local_op(0)->Rescale(0.0);
  op_acc_->AddAccumulationDelta(tcc_aux, wc_prev, wc, dt, "cell");

  UpdateBoundaryData(t_old, t_new, component);

  op_adv_->UpdateMatrices(S_->GetPtr<CompositeVector>(vol_flowrate_key_, Tags::DEFAULT).ptr());
  op_adv_->ApplyBCs(true, true, true);

  if (use_dispersion_) {
    int phase;
    double md;
    CalculateDispersionTensor_(*transport_phi, wc_c);
    FindDiffusionValue(component_names_[component], &md, &phase);
    if (md != 0.0) CalculateDiffusionTensor_(md, phase, *transport_phi, sat_c);

    op_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
    op_diff_->ApplyBCs(true, true, true);
  }

  // add sources
  Epetra_MultiVector& rhs_cell = *op_->rhs()->ViewComponent("cell");
  ComputeSources_(t_new, dt, rhs_cell, *tcc->ViewComponent("cell"), component, component);
}


/* ******************************************************************* 
* Boundary data
******************************************************************* */
void TransportImplicit_PK::UpdateBoundaryData(double t_old, double t_new, int component)
{
  for (int i = 0; i < bcs_.size(); i++) {
    bcs_[i]->Compute(t_old, t_new);
    bcs_[i]->ComputeSubmodel(mesh_, tcc);
  }

  auto& values = op_bc_->bc_value();
  auto& models = op_bc_->bc_model();

  ComputeBCs_(models, values, component);
}

}  // namespace Transport
}  // namespace Amazni
