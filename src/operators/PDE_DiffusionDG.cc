/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#include <vector>

// Amanzi
#include "DG_Modal.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"

// Operators
#include "InterfaceWhetStone.hh"
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
void
PDE_DiffusionDG::Init_(Teuchos::ParameterList& plist)
{
  Teuchos::ParameterList& schema_list = plist.sublist("schema");

  // parameters
  // -- discretization details
  Schema my_schema;
  auto base = my_schema.StringToKind(schema_list.get<std::string>("base"));

  matrix_ = plist.get<std::string>("matrix type");
  method_order_ = schema_list.get<int>("method order", 0);
  numi_order_ = plist.get<int>("quadrature order", method_order_);

  dg_ = Teuchos::rcp(new WhetStone::DG_Modal(schema_list, mesh_));
  my_schema.Init(dg_, mesh_, base);

  local_schema_col_ = my_schema;
  local_schema_row_ = my_schema;

  // create or check the existing Operator
  if (global_op_ == Teuchos::null) {
    global_schema_col_ = my_schema;
    global_schema_row_ = my_schema;

    // build the CVS from the global schema
    int nk;
    std::tie(std::ignore, std::ignore, nk) = *my_schema.begin();

    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true);
    cvs->AddComponent("cell", AmanziMesh::Entity_kind::CELL, nk);

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

  // -- continuity terms
  Schema tmp_schema;
  Teuchos::ParameterList schema_copy = schema_list;
  tmp_schema.Init(dg_, mesh_, AmanziMesh::Entity_kind::FACE);

  jump_up_op_ = Teuchos::rcp(new Op_Face_Schema(tmp_schema, tmp_schema, mesh_));
  global_op_->OpPushBack(jump_up_op_);

  jump_pu_op_ = Teuchos::rcp(new Op_Face_Schema(tmp_schema, tmp_schema, mesh_));
  global_op_->OpPushBack(jump_pu_op_);

  // -- stability jump term
  penalty_op_ = Teuchos::rcp(new Op_Face_Schema(tmp_schema, tmp_schema, mesh_));
  global_op_->OpPushBack(penalty_op_);
}


/* ******************************************************************
* Populate face-based 2x2 matrices on interior faces and 1x1 matrices
* on boundary faces.
****************************************************************** */
void
PDE_DiffusionDG::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                                const Teuchos::Ptr<const CompositeVector>& u)
{
  WhetStone::DenseMatrix Acell, Aface;

  double Kf(1.0);
  // volumetric term
  for (int c = 0; c != ncells_owned; ++c) {
    interface_->StiffnessMatrix(c, Acell);
    local_op_->matrices[c] = Acell;
  }

  // strenghen stability term
  for (int f = 0; f != nfaces_owned; ++f) {
    if (Kf_.get()) Kf = (*Kf_)[f];
    dg_->FaceMatrixPenalty(f, Kf, Aface);
    penalty_op_->matrices[f] = Aface;
  }

  // stability terms
  for (int f = 0; f != nfaces_owned; ++f) {
    auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    int c1 = cells[0];
    int c2 = (cells.size() > 1) ? cells[1] : c1;

    interface_->FaceMatrixJump(f, c1, c2, Aface);

    Aface *= -1.0;
    jump_up_op_->matrices[f] = Aface;

    Aface.Transpose();
    jump_pu_op_->matrices[f] = Aface;
  }
}


/* ******************************************************************
* Apply boundary conditions to the local matrices.
****************************************************************** */
void
PDE_DiffusionDG::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  AMANZI_ASSERT(bcs_trial_.size() > 0);

  const std::vector<int>& bc_model = bcs_trial_[0]->bc_model();
  const std::vector<std::vector<double>>& bc_value = bcs_trial_[0]->bc_value_vector();
  int nk = bc_value[0].size();

  Epetra_MultiVector& rhs_c = *global_op_->rhs()->ViewComponent("cell", true);

  int d = mesh_->getSpaceDimension();
  std::vector<AmanziGeometry::Point> tau(d - 1);

  // create integration object for all mesh cells
  WhetStone::NumericalIntegration numi(mesh_);

  for (int f = 0; f != nfaces_owned; ++f) {
    if (bc_model[f] == OPERATOR_BC_DIRICHLET || bc_model[f] == OPERATOR_BC_NEUMANN) {
      auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      int c = cells[0];

      const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);

      // set polynomial with Dirichlet data
      WhetStone::DenseVector coef(nk);
      for (int i = 0; i < nk; ++i) { coef(i) = bc_value[f][i]; }

      WhetStone::Polynomial pf(d, method_order_, coef);
      pf.set_origin(xf);

      // convert boundary polynomial to space polynomial
      pf.ChangeOrigin(mesh_->getCellCentroid(c));

      // extract coefficients and update right-hand side
      WhetStone::DenseMatrix& Pcell = penalty_op_->matrices[f];
      int nrows = Pcell.NumRows();
      int ncols = Pcell.NumCols();

      WhetStone::DenseVector v(nrows), pv(ncols), jv(ncols);
      v = pf.coefs();
      dg_->cell_basis(c).ChangeBasisNaturalToMy(v);

      if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
        WhetStone::DenseMatrix& Jcell = jump_up_op_->matrices[f];
        Pcell.Multiply(v, pv, false);
        Jcell.Multiply(v, jv, false);

        for (int i = 0; i < ncols; ++i) { rhs_c[i][c] += pv(i) + jv(i); }
      } else if (bc_model[f] == OPERATOR_BC_NEUMANN) {
        WhetStone::DenseMatrix& Jcell = jump_pu_op_->matrices[f];
        Jcell.Multiply(v, jv, false);

        for (int i = 0; i < ncols; ++i) { rhs_c[i][c] -= jv(i); }

        Pcell.PutScalar(0.0);
        Jcell.PutScalar(0.0);
        jump_up_op_->matrices[f].PutScalar(0.0);
      }
    }
  }
}


/* ******************************************************************
* Calculate mass flux from cell-centered data
****************************************************************** */
void
PDE_DiffusionDG::UpdateFlux(const Teuchos::Ptr<const CompositeVector>& solution,
                            const Teuchos::Ptr<CompositeVector>& flux)
{}

} // namespace Operators
} // namespace Amanzi
