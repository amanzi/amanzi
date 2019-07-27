/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

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
#include "Explicit_TI_RK.hh"
#include "FieldEvaluator.hh"
#include "GMVMesh.hh"
#include "LinearOperatorDefs.hh"
#include "LinearOperatorFactory.hh"
#include "Mesh.hh"
#include "PDE_Accumulation.hh"
#include "PDE_Advection.hh"
#include "PDE_AdvectionUpwind.hh"
#include "PDE_Diffusion.hh"
#include "PDE_DiffusionFactory.hh"
#include "PK_DomainFunctionFactory.hh"
#include "PK_Utils.hh"
#include "WhetStoneDefs.hh"

// amanzi::Transport
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
  std::string pk_name = pk_tree.name();
  auto found = pk_name.rfind("->");
  if (found != std::string::npos) pk_name.erase(0, found + 2);

  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(glist, "PKs", true);
  tp_list_ = Teuchos::sublist(pk_list, pk_name, true);
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
  vo_ = Teuchos::rcp(new VerboseObject("TransportImpl-" + domain, vo_list)); 

  // Create pointers to the primary flow field pressure.
  const auto& solution = S->GetFieldData(tcc_key_, "state");
  soln_->SetData(solution); 
  
  // boundary conditions
  op_bc_ = Teuchos::rcp(new Operators::BCs(mesh_, AmanziMesh::FACE, WhetStone::DOF_Type::SCALAR));

  auto values = op_bc_->bc_value();
  auto models = op_bc_->bc_model();

  for (int i = 0; i < bcs_.size(); i++) {
    std::vector<int>& tcc_index = bcs_[i]->tcc_index();
    int ncomp = tcc_index.size();
    
    for (auto bc = bcs_[i]->begin(); bc != bcs_[i]->end(); ++bc) {
      int f = bc->first;
      models[f] = Operators::OPERATOR_BC_DIRICHLET;
      std::vector<double>& bcval = bc->second;             
      values[f] = bcval[0];  // Only for one component at the moment
    }
  }

  Teuchos::ParameterList& oplist = tp_list_->sublist("operators")
                                            .sublist("advection operator")
                                            .sublist("matrix");

  op_adv_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(oplist, mesh_));
  op_adv_->SetBCs(op_bc_, op_bc_);
  
  op_ = op_adv_->global_operator();

  Teuchos::RCP<const CompositeVector> flux = S->GetFieldData(darcy_flux_key_);
  op_adv_->Setup(*flux);
  op_adv_->UpdateMatrices(flux.ptr());
  
  op_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op_));
  acc_term_ = Teuchos::rcp(new CompositeVector(*(S->GetFieldData(tcc_key_))));

  Epetra_MultiVector& acc_term_c = *acc_term_->ViewComponent("cell");

  for (int c = 0; c < ncells_owned; c++) {
    acc_term_c[0][c] = mesh_->cell_volume(c) * (*phi)[0][c];  // * (*ws_start)[0][c];
  }

  op_acc_->AddAccumulationTerm(*acc_term_, "cell");

  op_->SymbolicAssembleMatrix();
  op_->AssembleMatrix(); 

  // generic linear solver
  // AMANZI_ASSERT(ti_list_->isParameter("linear solver"));
  // solver_name_ = ti_list_->get<std::string>("linear solver");

  // preconditioner
  // AMANZI_ASSERT(ti_list_->isParameter("preconditioner"));
  // std::string name = ti_list_->get<std::string>("preconditioner");
  // Teuchos::ParameterList pc_list = preconditioner_list_->sublist(name);
  // op_->InitializePreconditioner(pc_list);

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << vo_->color("green") << "Initialization of PK is complete." 
               << vo_->reset() << std::endl << std::endl;
  }
}

}  // namespace Transport
}  // namespace Amazni
