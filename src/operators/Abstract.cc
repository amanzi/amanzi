/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

#include "MFD3D_Diffusion.hh"

#include "Abstract.hh"
#include "OperatorDefs.hh"
#include "Operator_Schema.hh"
#include "Op_Cell_Schema.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialize operator from parameter list.
****************************************************************** */
void Abstract::Init_(Teuchos::ParameterList& plist)
{
  Teuchos::ParameterList& range = plist.sublist("schema range");
  Teuchos::ParameterList& domain = plist.sublist("schema domain");

  if (global_op_ == Teuchos::null) {
    // constructor was given a mesh
    // -- range schema and cvs
    local_schema_row_.Init(range, mesh_);
    global_schema_row_ = local_schema_row_;

    Teuchos::RCP<CompositeVectorSpace> cvs_row = Teuchos::rcp(new CompositeVectorSpace());
    cvs_row->SetMesh(mesh_)->SetGhosted(true);

    for (auto it = global_schema_row_.items().begin(); it != global_schema_row_.items().end(); ++it) {
      std::string name(local_schema_row_.KindToString(it->kind));
      cvs_row->AddComponent(name, it->kind, it->num);
    }

    // -- domain schema and cvs
    local_schema_col_.Init(domain, mesh_);
    global_schema_col_ = local_schema_col_;

    Teuchos::RCP<CompositeVectorSpace> cvs_col = Teuchos::rcp(new CompositeVectorSpace());
    cvs_col->SetMesh(mesh_)->SetGhosted(true);

    for (auto it = global_schema_col_.items().begin(); it != global_schema_col_.items().end(); ++it) {
      std::string name(local_schema_col_.KindToString(it->kind));
      cvs_col->AddComponent(name, it->kind, it->num);
    }

    global_op_ = Teuchos::rcp(new Operator_Schema(cvs_row, cvs_col, plist, global_schema_row_, global_schema_col_));
    if (local_schema_col_.base() == AmanziMesh::CELL) {
      local_op_ = Teuchos::rcp(new Op_Cell_Schema(global_schema_row_, global_schema_col_, mesh_));
    }

  } else {
    // constructor was given an Operator
    global_schema_row_ = global_op_->schema_row();
    global_schema_col_ = global_op_->schema_col();

    mesh_ = global_op_->DomainMap().Mesh();
    local_schema_row_.Init(range, mesh_);
    local_schema_col_.Init(domain, mesh_);

    local_op_ = Teuchos::rcp(new Op_Cell_Schema(global_schema_row_, global_schema_col_, mesh_));
  }

  // register the advection Op
  global_op_->OpPushBack(local_op_);

  // parameters
  // -- discretization method
  method_ = plist.get<std::string>("discretization", "none");
}


/* ******************************************************************
* Loop over existing methods.
****************************************************************** */
void Abstract::UpdateMatrices()
{
  std::vector<WhetStone::DenseMatrix>& matrix = local_op_->matrices;
  std::vector<WhetStone::DenseMatrix>& matrix_shadow = local_op_->matrices_shadow;

  int dir, d(mesh_->space_dimension());
  AmanziMesh::Entity_ID_List nodes;

  WhetStone::MFD3D_Diffusion mfd(mesh_);
  WhetStone::DenseMatrix Wcell, Acell;

  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = 1.0;

  for (int c = 0; c < ncells_owned; ++c) {
    if (K_.get()) Kc = (*K_)[c];
    mfd.MassMatrixInverseGeneralized(c, Kc, Wcell);
    mfd.HybridizeGeneralized(c, Wcell, Acell);

    matrix[c] = Acell;
  }
}


/* *******************************************************************
* Apply boundary condition to the local matrices
******************************************************************* */
void Abstract::ApplyBCs(bool primary, bool eliminate)
{
  for (auto bc = bcs_trial_.begin(); bc != bcs_trial_.end(); ++bc) {
    if ((*bc)->kind() == AmanziMesh::FACE) {
      ApplyBCs_Face(bc->ptr(), local_op_, primary, eliminate);
    } 
    else if ((*bc)->kind() == AmanziMesh::NODE) {
      ApplyBCs_Node(bc->ptr(), local_op_, primary, eliminate);
    }
  }
}

}  // namespace Operators
}  // namespace Amanzi
