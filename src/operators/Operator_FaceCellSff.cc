/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  Operator whose unknowns are CELL + FACE, but which assembles the
  FACE only system and Schur complements cells.

  This uses special assembly.  Apply is done as if we had the full FACE+CELL
  system.  SymbolicAssembly() is done as if we had the CELL system, but with an
  additional step to get the layout due to the Schur'd system on FACE+CELL.
  Assemble, however, is done using a totally different approach.
*/

#include <vector>
#include "EpetraExt_RowMatrixOut.h"

#include "SuperMap.hh"
#include "GraphFE.hh"
#include "MatrixFE.hh"
#include "Preconditioner.hh"

#include "OperatorDefs.hh"
#include "OperatorUtils.hh"
#include "Op.hh"
#include "Op_Cell_FaceCell.hh"
#include "Op_Cell_Face.hh"
#include "Operator_FaceCellSff.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Special inverse algorithm is required to deal with Schur complement.
****************************************************************** */
int Operator_FaceCellSff::ApplyInverse(const CompositeVector& X, CompositeVector& Y) const
{
  // Check preconditions -- Sff must have exactly one CELL+FACE schema,
  // and no other CELL schemas that are not simply diagonal CELL_CELL.
  // Additionally, collect the diagonal for inversion.
  Epetra_MultiVector D_c(mesh_->cell_map(false),1);
  
  int num_with_cells = 0;
  for (const_op_iterator it = begin(); it != end(); ++it) {
    if ((*it)->schema_old_ & OPERATOR_SCHEMA_DOFS_CELL) {
      if (((*it)->schema_old_ == (OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_CELL))
          && ((*it)->diag->MyLength() == ncells_owned)) {
        // diagonal schema
        for (int c = 0; c != ncells_owned; ++c) {
          D_c[0][c] += (*(*it)->diag)[0][c];
        }
      } else {
        num_with_cells++;
      }
    }
  }

  // This error is a false positive when doing coupled surface/subsurface runs. --etc
  // if (num_with_cells == 0 || num_with_cells > 1) {
  //   Errors::Message msg("Schur complement to Sff must have exactly one cell-based schema with off-diagonal entries or other, non-cell dofs.");
  //   Exceptions::amanzi_throw(msg);
  // }

  int ierr(0);
  Y.PutScalarGhosted(0.0);

  // apply preconditioner inversion
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell");
  // const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

  // Temporary cell and face vectors.
  CompositeVector T(X);
  T.PutScalarGhosted(0.0);

  for (const_op_iterator it = begin(); it != end(); ++it) {
    if ((*it)->schema_old_ == (OPERATOR_SCHEMA_BASE_CELL |
                               OPERATOR_SCHEMA_DOFS_CELL | OPERATOR_SCHEMA_DOFS_FACE)) {

      // FORWARD ELIMINATION:  Tf = Xf - Afc inv(Acc) Xc
      AmanziMesh::Entity_ID_List faces;
      {
        Epetra_MultiVector& Tf = *T.ViewComponent("face", true);
        for (int c = 0; c < ncells_owned; c++) {
          mesh_->cell_get_faces(c, &faces);
          int nfaces = faces.size();

          WhetStone::DenseMatrix& Acell = (*it)->matrices[c];

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
        const Epetra_MultiVector& Tf_short = *T.ViewComponent("face", false);
        Epetra_MultiVector& Yf_short = *Y.ViewComponent("face", false);

        ierr = preconditioner_->ApplyInverse(Tf_short, Yf_short);
        AMANZI_ASSERT(ierr >= 0);
      }

      Y.ScatterMasterToGhosted("face");

      {
        const Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);
        Epetra_MultiVector& Yc = *Y.ViewComponent("cell", false);
        // BACKWARD SUBSTITUTION:  Yc = inv(Acc) (Xc - Acf Yf)
        for (int c = 0; c < ncells_owned; c++) {
          mesh_->cell_get_faces(c, &faces);
          int nfaces = faces.size();

          WhetStone::DenseMatrix& Acell = (*it)->matrices[c];

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


/* ******************************************************************
* Special assemble algorithm is required to deal with Schur complement.
****************************************************************** */
void Operator_FaceCellSff::AssembleMatrix(const SuperMap& map, MatrixFE& matrix,
                                          int my_block_row, int my_block_col) const
{
  // Check preconditions -- Scc must have exactly one CELL+FACE schema,
  // and no other CELL schemas that are not simply diagonal CELL_CELL.
  // Additionally, collect the diagonal for inversion.
  Epetra_MultiVector D_c(mesh_->cell_map(false), 1);
  
  int num_with_cells = 0;
  for (const_op_iterator it = begin(); it != end(); ++it) {
    if ((*it)->schema_old_ & OPERATOR_SCHEMA_DOFS_CELL) {
      if (((*it)->schema_old_ == (OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_CELL))
          && ((*it)->diag->MyLength() == ncells_owned)) {
        // diagonal schema
        for (int c = 0; c != ncells_owned; ++c) {
          D_c[0][c] += (*(*it)->diag)[0][c];
        }
      } else {
        num_with_cells++;
      }
    }
  }

  // This error is a false positive when doing coupled surface/subsurface runs. --etc
  // if (num_with_cells == 0 || num_with_cells > 1) {
  //   Errors::Message msg("Schur complement to Sff must have exactly one cell-based schema with off-diagonal entries or other, non-cell dofs.");
  //   Exceptions::amanzi_throw(msg);
  // }
  
  // schur complement
  int i_schur = 0;
  for (const_op_iterator it = begin(); it != end(); ++it) {
    if ((*it)->schema_old_ == (OPERATOR_SCHEMA_BASE_CELL |
                               OPERATOR_SCHEMA_DOFS_CELL | OPERATOR_SCHEMA_DOFS_FACE)) {
      AMANZI_ASSERT((*it)->matrices.size() == ncells_owned);

      // create or get extra ops, and keep them for future use
      Teuchos::RCP<Op_Cell_Face> schur_op;

      if (schur_ops_.size() > i_schur) {
        // get existing
        schur_op = schur_ops_[i_schur];
      } else {
        // create and fill
        AMANZI_ASSERT(i_schur == schur_ops_.size());
        std::string name = "Sff alt CELL_FACE";
        schur_op = Teuchos::rcp(new Op_Cell_Face(name, mesh_));
        schur_ops_.push_back(schur_op);
        for (int c = 0; c != ncells_owned; ++c) {
          int nfaces = mesh_->cell_get_num_faces(c);
          schur_op->matrices[c] = WhetStone::DenseMatrix(nfaces, nfaces);
        }
      }

      // populate the schur component
      AmanziMesh::Entity_ID_List faces;
      for (int c = 0; c != ncells_owned; ++c) {
        mesh_->cell_get_faces(c, &faces);
        int nfaces = faces.size();

        WhetStone::DenseMatrix& Scell = schur_op->matrices[c];
        WhetStone::DenseMatrix& Acell = (*it)->matrices[c];

        double tmp = Acell(nfaces, nfaces) + D_c[0][c];

        for (int n = 0; n < nfaces; n++) {
          for (int m = 0; m < nfaces; m++) {
            Scell(n, m) = Acell(n, m) - Acell(n, nfaces) * Acell(nfaces, m) / tmp;
          }
        }
      }

      // Assemble this Schur Op into matrix
      schur_op->AssembleMatrixOp(this, map, matrix, my_block_row, my_block_col);
    } else if (((*it)->schema_old_ ==
                (OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_CELL))
               && ((*it)->diag->MyLength() == ncells_owned)) {
      // pass, already part of cell inverse
    } else {
      (*it)->AssembleMatrixOp(this, map, matrix, my_block_row, my_block_col);
    }
  }
}


/* ******************************************************************
* Visit method for Apply -- this is identical to Operator_FaceCell's
* version.
****************************************************************** */
int Operator_FaceCellSff::ApplyMatrixFreeOp(const Op_Cell_FaceCell& op,
                                            const CompositeVector& X, CompositeVector& Y) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

  X.ScatterMasterToGhosted();
  const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell");
  Y.PutScalarGhosted(0.);
  
  {
    Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);
    Epetra_MultiVector& Yc = *Y.ViewComponent("cell");

    AmanziMesh::Entity_ID_List faces;
    for (int c=0; c!=ncells_owned; ++c) {
      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      WhetStone::DenseVector v(nfaces + 1), av(nfaces + 1);
      for (int n=0; n!=nfaces; ++n) {
        v(n) = Xf[0][faces[n]];
      }
      v(nfaces) = Xc[0][c];

      const WhetStone::DenseMatrix& Acell = op.matrices[c];
      Acell.Multiply(v, av, false);

      for (int n=0; n!=nfaces; ++n) {
        Yf[0][faces[n]] += av(n);
      }
      Yc[0][c] += av(nfaces);
    }
  }

  Y.GatherGhostedToMaster(Add);
  return 0;
}


/* ******************************************************************
* Create a global matrix.
****************************************************************** */
void Operator_FaceCellSff::SymbolicAssembleMatrix()
{
  // SuperMap for Sff is face only
  CompositeVectorSpace smap_space;
  smap_space.SetMesh(mesh_)->SetComponent("face", AmanziMesh::FACE, 1);
  smap_ = createSuperMap(smap_space);

  // create the graph
  int row_size = MaxRowSize(*mesh_, schema(), 1);
  Teuchos::RCP<GraphFE> graph = Teuchos::rcp(new GraphFE(smap_->Map(),
          smap_->GhostedMap(), smap_->GhostedMap(), row_size));

  Operator::SymbolicAssembleMatrix(*smap_, *graph, 0, 0);

  // Completing and optimizing the graphs
  int ierr = graph->FillComplete(smap_->Map(), smap_->Map());
  AMANZI_ASSERT(!ierr);

  // create global matrix
  Amat_ = Teuchos::rcp(new MatrixFE(graph));
  A_ = Amat_->Matrix();
}


/* ******************************************************************
* visit method for sparsity structure of Schur complement
****************************************************************** */
void Operator_FaceCellSff::SymbolicAssembleMatrixOp(const Op_Cell_FaceCell& op,
                                                    const SuperMap& map, GraphFE& graph,
                                                    int my_block_row, int my_block_col) const
{
  std::string name = "Sff alt CELL_FACE";
  Op_Cell_Face schur_op(name, mesh_);
  Operator_FaceCell::SymbolicAssembleMatrixOp(schur_op, map, graph, my_block_row, my_block_col);
}

}  // namespace Operators
}  // namespace Amanzi

