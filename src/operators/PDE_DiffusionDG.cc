/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

// Amanzi
#include "MFD3DFactory.hh"

// Operators
#include "Op.hh"
#include "Op_Cell_Schema.hh"
#include "Op_Face_Schema.hh"
#include "OperatorDefs.hh"
#include "Operator_Schema.hh"
#include "PDE_DiffusionDG.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization
****************************************************************** */
void PDE_DiffusionDG::Init_(Teuchos::ParameterList& plist)
{
  // create the local Op and register it with the global Operator
  Schema my_schema;
  Teuchos::ParameterList& schema_list = plist.sublist("schema");
  my_schema.Init(schema_list, mesh_);

  local_schema_col_ = my_schema;
  local_schema_row_ = my_schema;

  // create or check the existing Operator
  if (global_op_ == Teuchos::null) {
    global_schema_col_ = my_schema;
    global_schema_row_ = my_schema;

    // build the CVS from the global schema
    int nk = my_schema.begin()->num;
    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true);
    cvs->AddComponent("cell", AmanziMesh::CELL, nk);

    global_op_ = Teuchos::rcp(new Operator_Schema(cvs, plist, my_schema));

  } else {
    // constructor was given an Operator
    global_schema_col_ = global_op_->schema_col();
    global_schema_row_ = global_op_->schema_row();
    mesh_ = global_op_->DomainMap().Mesh();
  }

  // create local operators
  // -- stiffness
  local_op_ = Teuchos::rcp(new Op_Cell_Schema(my_schema, my_schema, mesh_));
  global_op_->OpPushBack(local_op_);

  // -- jump
  Teuchos::ParameterList schema_copy = schema_list;
  schema_copy.set<std::string>("base", "face");
  my_schema.Init(schema_copy, mesh_);

  jump_uu_op_ = Teuchos::rcp(new Op_Face_Schema(my_schema, my_schema, mesh_));
  global_op_->OpPushBack(jump_uu_op_);

  // parameters
  // -- discretization details
  method_ = plist.get<std::string>("method");
  method_order_ = plist.get<int>("method order", 0);
  matrix_ = plist.get<std::string>("matrix type");
}


/* ******************************************************************
* Setup methods: scalar coefficients
****************************************************************** */
void PDE_DiffusionDG::SetProblemCoefficients(
   const std::shared_ptr<std::vector<WhetStone::Tensor> >& Kc,
   const std::shared_ptr<std::vector<double> >& Kf) 
{
  Kc_ = Kc;
  Kf_ = Kf;
}


/* ******************************************************************
* Populate face-based 2x2 matrices on interior faces and 1x1 matrices
* on boundary faces.
****************************************************************** */
void PDE_DiffusionDG::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                     const Teuchos::Ptr<const CompositeVector>& u)
{
  WhetStone::DG_Modal dg(method_order_, mesh_);
  WhetStone::DenseMatrix Acell, Aface;

  double Kf(1.0);
  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = 1.0;

  for (int c = 0; c != ncells_owned; ++c) {
    if (Kc_.get()) Kc = (*Kc_)[c];
    dg.StiffnessMatrix(c, Kc, Acell);
    local_op_->matrices[c] = Acell;
  }

  for (int f = 0; f != nfaces_owned; ++f) {
    if (Kf_.get()) Kf = (*Kf_)[f];
    dg.JumpMatrix(f, Kf, Aface);
    jump_uu_op_->matrices[f] = Aface;
  }
}


/* ******************************************************************
* Apply boundary conditions to the local matrices.
****************************************************************** */
void PDE_DiffusionDG::ApplyBCs(bool primary, bool eliminate)
{
  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  const std::vector<std::vector<double> >& bc_value = bcs_trial_[0]->bc_value_vector();
  int d = bc_value[0].size();

  AmanziMesh::Entity_ID_List cells;

  Epetra_MultiVector& rhs_c = *global_op_->rhs()->ViewComponent("cell", true);
  const Schema& schema = global_op_->schema_row();

  for (int f = 0; f != nfaces_owned; ++f) {
    if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
      mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
      int c = cells[0];

      WhetStone::DenseMatrix& Acell = jump_uu_op_->matrices[f];
      int nrows = Acell.NumRows();
      int ncols = Acell.NumCols();

      const std::vector<double>& value = bc_value[f];
      WhetStone::DenseVector v(nrows), av(ncols);

      v.PutScalar(0.0);
      Acell.Multiply(v, av, false);

      for (int i = 0; i < ncols; ++i) {
        rhs_c[i][c] += av(i);
      }
    }
  } 
}


/* ******************************************************************
* Calculate mass flux from cell-centered data
****************************************************************** */
void PDE_DiffusionDG::UpdateFlux(const Teuchos::Ptr<const CompositeVector>& solution,
                                 const Teuchos::Ptr<CompositeVector>& flux)
{
}

}  // namespace Operators
}  // namespace Amanzi


