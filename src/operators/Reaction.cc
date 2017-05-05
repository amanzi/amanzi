/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "dg.hh"

#include "Op_Cell_Schema.hh"
#include "OperatorDefs.hh"
#include "Operator_Schema.hh"
#include "Reaction.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialize operator from parameter list.
****************************************************************** */
void Reaction::InitReaction_(Teuchos::ParameterList& plist)
{
  Teuchos::ParameterList& range = plist.sublist("schema");

  if (global_op_ == Teuchos::null) {
    // constructor was given a mesh
    local_schema_row_.Init(range, mesh_);
    global_schema_row_ = local_schema_row_;

    local_schema_col_ = local_schema_row_;
    global_schema_col_ = global_schema_row_;

    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true);

    for (auto it = global_schema_row_.items().begin(); it != global_schema_row_.items().end(); ++it) {
      std::string name(local_schema_row_.KindToString(it->kind));
      cvs->AddComponent(name, it->kind, it->num);
    }

    global_op_ = Teuchos::rcp(new Operator_Schema(cvs, cvs, plist, global_schema_row_, global_schema_col_));
    local_op_ = Teuchos::rcp(new Op_Cell_Schema(global_schema_row_, global_schema_col_, mesh_));

  } else {
    // constructor was given an Operator
    global_schema_row_ = global_op_->schema_row();
    global_schema_col_ = global_op_->schema_col();

    mesh_ = global_op_->DomainMap().Mesh();
    local_schema_row_.Init(range, mesh_);
    local_schema_col_ = local_schema_row_;

    local_op_ = Teuchos::rcp(new Op_Cell_Schema(global_schema_row_, global_schema_col_, mesh_));
  }

  // register the advection Op
  global_op_->OpPushBack(local_op_);

  // parameters
  std::string name = plist.get<std::string>("discretization", "none");
  if (name == "DG order 0") {
    space_col_ = DG0;
    space_row_ = DG0;
  } else if (name == "DG order 1") {
    space_col_ = DG1;
    space_row_ = DG1;
  } else {
    Errors::Message msg;
    msg << "Discretization method is either missing or invalid.";
    Exceptions::amanzi_throw(msg);
  }

  // mesh info
  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);

  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::USED);
}


/* ******************************************************************
* Collection of local matrices.
****************************************************************** */
void Reaction::UpdateMatrices(const CompositeVector& u)
{
  std::vector<WhetStone::DenseMatrix>& matrix = local_op_->matrices;
  std::vector<WhetStone::DenseMatrix>& matrix_shadow = local_op_->matrices_shadow;

  AmanziMesh::Entity_ID_List nodes;
  int d = mesh_->space_dimension();

  WhetStone::DG dg(mesh_);

  int k(0);
  if (space_col_ == DG0 && space_row_ == DG0) { 
    k = 0;
  } else if (space_col_ == DG1 && space_row_ == DG1) {  
    k = 1;
  }

  for (int c = 0; c < ncells_owned; ++c) {
    int ndofs = (k + 2) * (k + 1) / 2;

    WhetStone::DenseMatrix Mcell(ndofs, ndofs);
    dg.TaylorMassMatrix(c, k, Mcell);

    matrix[c] = Mcell;
  }
}


/* *******************************************************************
* Apply boundary condition to the local matrices
******************************************************************* */
void Reaction::ApplyBCs(bool primary, bool eliminate)
{
  for (auto bc = bcs_trial_.begin(); bc != bcs_trial_.end(); ++bc) {
  }
}


/* *******************************************************************
* Identify the advected flux of u
******************************************************************* */
void Reaction::UpdateFlux(
    const CompositeVector& h, const CompositeVector& u,
    const Teuchos::RCP<BCs>& bc, CompositeVector& flux)
{
}

}  // namespace Operators
}  // namespace Amanzi
