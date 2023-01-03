/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! Provides ApplyInverse() using a Schur complement.
#include "MeshDefs.hh"
#include "Op.hh"
#include "InverseFactory.hh"
#include "InverseSchurComplement.hh"

namespace Amanzi {
namespace AmanziSolvers {

void
InverseSchurComplement::set_inverse_parameters(Teuchos::ParameterList& plist)
{
  std::string method_name = plist.get<std::string>("preconditioning method");
  solver_ = createAssembledMethod<>(method_name, plist);
}

void
InverseSchurComplement::InitializeInverse()
{
  AMANZI_ASSERT(h_.get());      // set_matrices was called
  AMANZI_ASSERT(solver_.get()); // set_inverse_parameters was called
  solver_->set_matrix(Teuchos::null);
  h_->SymbolicAssembleMatrix();
  solver_->set_matrix(h_->A());
  solver_->InitializeInverse();
}

void
InverseSchurComplement::ComputeInverse()
{
  h_->AssembleMatrix();
  solver_->ComputeInverse();
}

int
InverseSchurComplement::ApplyInverse(const CompositeVector& X, CompositeVector& Y) const
{
  // Check preconditions -- Sff must have exactly one "coupling" CELL+FACE
  // schema which connects cells to faces, and no other CELL schemas that are
  // not simply diagonal CELL_CELL.  Additionally, collect the diagonal for
  // inversion.
  Epetra_MultiVector D_c(h_->Mesh()->getMap(AmanziMesh::Entity_kind::CELL,false), 1);
  int ncells_owned = D_c.MyLength();
  int nfaces_owned = h_->Mesh()->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_type::OWNED);

  int nschema_coupling = 0;
  for (const auto& op : *h_) {
    if (op->Matches(Operators::OPERATOR_SCHEMA_DOFS_CELL, Operators::OPERATOR_SCHEMA_RULE_SUBSET)) {
      // HAS CELLS
      if (op->Matches(Operators::OPERATOR_SCHEMA_BASE_CELL | Operators::OPERATOR_SCHEMA_DOFS_CELL,
                      Operators::OPERATOR_SCHEMA_RULE_EXACT) &&
          (op->diag->MyLength() == ncells_owned)) {
        // diagonal schema
        for (int c = 0; c != ncells_owned; ++c) { D_c[0][c] += (*op->diag)[0][c]; }
      } else if (op->Matches(Operators::OPERATOR_SCHEMA_BASE_FACE,
                             Operators::OPERATOR_SCHEMA_RULE_SUBSET) &&
                 op->matrices.size() == nfaces_owned) {
        Errors::Message msg;
        msg << "Schema " << op->schema_string << " cannot be Schur complemented.";
        Exceptions::amanzi_throw(msg);
      } else if (op->Matches(Operators::OPERATOR_SCHEMA_BASE_CELL |
                               Operators::OPERATOR_SCHEMA_DOFS_CELL |
                               Operators::OPERATOR_SCHEMA_DOFS_FACE,
                             Operators::OPERATOR_SCHEMA_RULE_EXACT) &&
                 (op->matrices.size() == ncells_owned)) {
        nschema_coupling++;
      }
    }
  }
  if (nschema_coupling == 0) {
    Errors::Message msg("Schur complement has no FACE+CELL schema.");
    Exceptions::amanzi_throw(msg);
  } else if (nschema_coupling > 1) {
    Errors::Message msg("Schur complement has multiple FACE+CELL schema.");
    Exceptions::amanzi_throw(msg);
  }

  int ierr(0);
  Y.PutScalarGhosted(0.0);

  // apply preconditioner inversion
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell");

  // Temporary cell and face vectors.
  CompositeVector T(X);
  T.PutScalarGhosted(0.0);

  for (const auto& op : *h_) {
    if (op->Matches(Operators::OPERATOR_SCHEMA_BASE_CELL | Operators::OPERATOR_SCHEMA_DOFS_CELL |
                      Operators::OPERATOR_SCHEMA_DOFS_FACE,
                    Operators::OPERATOR_SCHEMA_RULE_EXACT) &&
        op->matrices.size() == ncells_owned) {
      // FORWARD ELIMINATION:  Tf = Xf - Afc inv(Acc) Xc
      {
        Epetra_MultiVector& Tf = *T.ViewComponent("face", true);
        for (int c = 0; c < ncells_owned; c++) {
          const auto& faces = h_->Mesh()->getCellFaces(c);
          int nfaces = faces.size();

          WhetStone::DenseMatrix& Acell = op->matrices[c];

          double tmp = Xc[0][c] / (Acell(nfaces, nfaces) + D_c[0][c]);
          for (int n = 0; n < nfaces; n++) {
            int f = faces[n];
            Tf[0][f] -= Acell(n, nfaces) * tmp;
          }
        }
      }

      T.GatherGhostedToMaster("face", Add);

      // Solve the Schur complement system Sff * Yf = Tf.
      {
        Epetra_MultiVector& Yf = *Y.ViewComponent("face", false);
        const Epetra_MultiVector& Tf = *T.ViewComponent("face", false);
        ierr = solver_->ApplyInverse(*Tf(0), *Yf(0));
        AMANZI_ASSERT(ierr >= 0);
      }

      Y.ScatterMasterToGhosted("face");

      {
        const Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);
        Epetra_MultiVector& Yc = *Y.ViewComponent("cell", false);
        // BACKWARD SUBSTITUTION:  Yc = inv(Acc) (Xc - Acf Yf)
        for (int c = 0; c < ncells_owned; c++) {
          const auto& faces = h_->Mesh()->getCellFaces(c);
          int nfaces = faces.size();

          WhetStone::DenseMatrix& Acell = op->matrices[c];

          double tmp = Xc[0][c];
          for (int n = 0; n < nfaces; n++) {
            int f = faces[n];
            tmp -= Acell(nfaces, n) * Yf[0][f];
          }
          Yc[0][c] = tmp / (Acell(nfaces, nfaces) + D_c[0][c]);
        }
      }
    }
  }
  return ierr;
}


} // namespace AmanziSolvers
} // namespace Amanzi
