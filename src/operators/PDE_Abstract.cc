/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include <vector>

#include "errors.hh"
#include "BilinearFormFactory.hh"

#include "Operator_Schema.hh"
#include "Op_Cell_Schema.hh"
#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
 * Initialize operator from parameter list.
 ****************************************************************** */
void
PDE_Abstract::Init_(Teuchos::ParameterList& plist)
{
  Teuchos::ParameterList range, domain;
  if (plist.isSublist("schema range") && plist.isSublist("schema domain")) {
    range = plist.sublist("schema range");
    domain = plist.sublist("schema domain");
  } else if (plist.isSublist("schema")) {
    range = plist.sublist("schema");
    domain = range;
  } else {
    Errors::Message msg;
    msg << "Schema mismatch for abstract operator.\n"
        << "  Use \"schema\" for a square operator.\n"
        << "  Use \"schema range\" and \"schema domain\" for a general "
           "operator.\n";
    Exceptions::amanzi_throw(msg);
  }

  if (global_op_ == Teuchos::null) {
    // constructor was given a mesh
    // -- range schema and cvs
    local_schema_row_.Init(range, mesh_);
    global_schema_row_ = local_schema_row_;

    Teuchos::RCP<CompositeVectorSpace> cvs_row =
      Teuchos::rcp(new CompositeVectorSpace());
    cvs_row->SetMesh(mesh_)->SetGhosted(true);

    for (auto it = global_schema_row_.begin(); it != global_schema_row_.end();
         ++it) {
      std::string name(local_schema_row_.KindToString(it->kind));
      cvs_row->AddComponent(name, it->kind, it->num);
    }

    // -- domain schema and cvs
    local_schema_col_.Init(domain, mesh_);
    global_schema_col_ = local_schema_col_;

    Teuchos::RCP<CompositeVectorSpace> cvs_col =
      Teuchos::rcp(new CompositeVectorSpace());
    cvs_col->SetMesh(mesh_)->SetGhosted(true);

    for (auto it = global_schema_col_.begin(); it != global_schema_col_.end();
         ++it) {
      std::string name(local_schema_col_.KindToString(it->kind));
      cvs_col->AddComponent(name, it->kind, it->num);
    }

    global_op_ = Teuchos::rcp(new Operator_Schema(
      cvs_row, cvs_col, plist, global_schema_row_, global_schema_col_));
    if (local_schema_col_.base() == AmanziMesh::CELL) {
      local_op_ = Teuchos::rcp(
        new Op_Cell_Schema(global_schema_row_, global_schema_col_, mesh_));
    }

  } else {
    // constructor was given an Operator
    global_schema_row_ = global_op_->schema_row();
    global_schema_col_ = global_op_->schema_col();

    mesh_ = global_op_->DomainMap().Mesh();
    local_schema_row_.Init(range, mesh_);
    local_schema_col_.Init(domain, mesh_);

    local_op_ = Teuchos::rcp(
      new Op_Cell_Schema(global_schema_row_, global_schema_col_, mesh_));
  }

  // register the advection Op
  global_op_->OpPushBack(local_op_);

  // parse parameters
  // -- discretization details
  matrix_ = plist.get<std::string>("matrix type");
  grad_on_test_ = plist.get<bool>("gradient operator on test function", true);

  // discretization method
  mfd_ = WhetStone::BilinearFormFactory::Create(plist, mesh_);
}


/* ******************************************************************
 * Populate containers of elemental matrices using MFD factory.
 * NOTE: input parameters are not yet used.
 ****************************************************************** */
void
PDE_Abstract::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                             const Teuchos::Ptr<const CompositeVector>& p)
{
  std::vector<WhetStone::DenseMatrix>& matrix = local_op_->matrices;
  std::vector<WhetStone::DenseMatrix>& matrix_shadow =
    local_op_->matrices_shadow;

  int dir, d(mesh_->space_dimension());
  AmanziMesh::Entity_ID_List nodes;

  // identify type of coefficient
  std::string coef("constant");
  if (Kpoly_.get()) coef = "polynomial";
  if (Kvec_.get()) coef = "vector polynomial";

  WhetStone::DenseMatrix Mcell, Acell;
  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = 1.0;

  if (matrix_ == "mass" && coef == "constant") {
    for (int c = 0; c < ncells_owned; ++c) {
      if (K_.get()) Kc = (*K_)[c];
      mfd_->MassMatrix(c, Kc, Acell);
      matrix[c] = Acell;
    }
  } else if (matrix_ == "mass" && coef == "polynomial") {
    for (int c = 0; c < ncells_owned; ++c) {
      mfd_->MassMatrix(c, (*Kpoly_)[c], Acell);
      matrix[c] = Acell;
    }
  } else if (matrix_ == "mass inverse" && coef == "constant") {
    for (int c = 0; c < ncells_owned; ++c) {
      if (K_.get()) Kc = (*K_)[c];
      mfd_->MassMatrixInverse(c, Kc, Acell);
      matrix[c] = Acell;
    }
  } else if (matrix_ == "mass inverse" && coef == "polynomial") {
    for (int c = 0; c < ncells_owned; ++c) {
      mfd_->MassMatrixInverse(c, (*Kpoly_)[c], Acell);
      matrix[c] = Acell;
    }
  } else if (matrix_ == "stiffness" && coef == "constant") {
    for (int c = 0; c < ncells_owned; ++c) {
      if (K_.get()) Kc = (*K_)[c];
      mfd_->StiffnessMatrix(c, Kc, Acell);
      matrix[c] = Acell;
    }
  } else if (matrix_ == "divergence") {
    for (int c = 0; c < ncells_owned; ++c) {
      mfd_->DivergenceMatrix(c, Acell);
      matrix[c] = Acell;
    }
  } else if (matrix_ == "advection" && coef == "constant") {
    const Epetra_MultiVector& u_c = *u->ViewComponent("cell", false);

    for (int c = 0; c < ncells_owned; ++c) {
      AmanziGeometry::Point vc(d);
      for (int i = 0; i < d; ++i) vc[i] = u_c[i][c];

      mfd_->AdvectionMatrix(c, vc, Acell, grad_on_test_);
      matrix[c] = Acell;
    }
  } else if (matrix_ == "advection" && coef == "vector polynomial") {
    for (int c = 0; c < ncells_owned; ++c) {
      mfd_->AdvectionMatrix(c, (*Kvec_)[c], Acell, grad_on_test_);
      matrix[c] = Acell;
    }
  } else {
    Errors::Message msg;
    msg << "Unsupported combination matrix=" << matrix_ << " and coef=" << coef
        << "\n";
    Exceptions::amanzi_throw(msg);
  }
}

} // namespace Operators
} // namespace Amanzi
