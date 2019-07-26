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
}


void TransportImplicit_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  Transport_PK::Initialize(S);

  Teuchos::ParameterList& oplist = tp_list_->sublist("operators")
                                            .sublist("advection operator")
                                            .sublist("matrix");

  op_adv_ = Teuchos::rcp(new Operators::PDE_AdvectionUpwind(oplist, mesh_));
  op_ = op_adv_->global_operator();

  Teuchos::RCP<const CompositeVector> flux = S->GetFieldData(darcy_flux_key_);
  op_adv_->Setup(*flux);
  op_adv_->UpdateMatrices(flux.ptr());
  
  op_acc_ = Teuchos::rcp(new Operators::PDE_Accumulation(AmanziMesh::CELL, op_));
  acc_term_ = Teuchos::rcp(new CompositeVector(*(S->GetFieldData(tcc_key_))));

  Epetra_MultiVector& acc_term_c = *acc_term_ -> ViewComponent("cell");

  for (int c=0; c<ncells_owned; c++){
    acc_term_c[0][c] = mesh_->cell_volume(c) * (*phi)[0][c] * (*ws_start)[0][c];
  }

  op_acc_->AddAccumulationTerm(*acc_term_, "cell");

  op_->SymbolicAssembleMatrix();
  op_->CreateCheckPoint();

  op_->AssembleMatrix(); 

  // // -- generic linear solver.
  // AMANZI_ASSERT(ti_list_->isParameter("linear solver"));
  // solver_name_ = ti_list_->get<std::string>("linear solver");

  // // -- preconditioner. There is no need to enhance it for Darcy
  // AMANZI_ASSERT(ti_list_->isParameter("preconditioner"));
  // std::string name = ti_list_->get<std::string>("preconditioner");
  // Teuchos::ParameterList pc_list = preconditioner_list_->sublist(name);
  // op_->InitializePreconditioner(pc_list);
}

}
}
