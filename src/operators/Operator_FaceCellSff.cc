/*
  This is the Operator component of the Amanzi code.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
  Ethan Coon (ecoon@lanl.gov)
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

/* ******************************************************************
   Operator whose unknowns are CELL + FACE, but which assembles the CELL only
   system and Schur complements the face.

   This uses special assembly.  Apply is done as if we had the full FACE+CELL
   system.  SymbolicAssembly() is done as if we had the CELL system, but with an
   additional step to get the layout due to the Schur'd system on FACE+CELL.
   Assemble, however, is done using a totally different approach.

   ---------------------------------------------------------------------

   1. Operator is a linear operator acting from linear space X to linear
   space Y. These spaces are described by CompositeVectors (CV). A few
   maps X->Y is supported.

   At the moment X = Y. Extension to TreeVectors should not be done in
   this class.

   2. Operator is an un-ordered additive collection of lower-rank (or
   equal) simple operators. During its construction, an operator can
   only grow by assimilating more operators.

   At the moment, an operator cannot be split into two operators, but
   there are no desing restriction for doing it in the future.

   3. A simple operator (a set of 1 operators) is defined by triple:
   scheme + elemental matrices + diagonal matrix. The schema specifies
   structure of elemental matrices, e.g. cell-based matrices
   representing interation between face-based unknowns.

   4. Operator can be converted to Epetra_FECrsMatrix matrix to generate
   a preconditioner. This operation cannot be applied to a subset of
   defining operators.

   Note. The operators can be initialized from other operators.
   Since data are never copied by default, we have to track
   down the ownership of data.
   ****************************************************************** */

namespace Amanzi {
namespace Operators {

// Special Apply Inverse required to deal with schur complement
int
Operator_FaceCellSff::ApplyInverse(const CompositeVector& X, CompositeVector& Y) const
{
  // Check preconditions -- Sff must have exactly one CELL+FACE schema,
  // and no other CELL schemas that are not simply diagonal CELL_CELL.
  // Additionally, collect the diagonal for inversion.
  Epetra_MultiVector D_c(mesh_->cell_map(false),1);
  
  int num_with_cells = 0;
  for (const_op_iterator it = OpBegin(); it != OpEnd(); ++it) {
    if ((*it)->schema & OPERATOR_SCHEMA_DOFS_CELL) {
      if ((*it)->schema == (OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_CELL)) {
        // diagonal schema
        for (int c = 0; c != ncells_owned; ++c) {
          D_c[0][c] += (*it)->vals[c];
        }
      } else {
        num_with_cells++;
      }
    }
  }
  if (num_with_cells == 0 || num_with_cells > 1) {
    Errors::Message msg("Schur complement to Sff must have exactly one cell-based schema with off-diagonal entries or other, non-cell dofs.");
    Exceptions::amanzi_throw(msg);
  }


  int ierr(0);
  Y.PutScalarGhosted(0.);

  // apply preconditioner inversion
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell");
  const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

  // Temporary cell and face vectors.
  CompositeVector T(X);
  T.PutScalarGhosted(0.);

  for (const_op_iterator it = OpBegin(); it != OpEnd(); ++it) {
    if ((*it)->schema == (OPERATOR_SCHEMA_BASE_CELL |
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

        // // dump the schur complement
        // std::stringstream filename_s2;
        // filename_s2 << "schur_PC_" << 0 << ".txt";
        // EpetraExt::RowMatrixToMatlabFile(filename_s2.str().c_str(), *A_);

        ierr = preconditioner_->ApplyInverse(Tf_short, Yf_short);
        ASSERT(!ierr);
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

// Special AssembleMatrix required to deal with schur complement
void
Operator_FaceCellSff::AssembleMatrix(const SuperMap& map, MatrixFE& matrix,
        int my_block_row, int my_block_col) const
{
  // Check preconditions -- Scc must have exactly one CELL+FACE schema,
  // and no other CELL schemas that are not simply diagonal CELL_CELL.
  // Additionally, collect the diagonal for inversion.
  Epetra_MultiVector D_c(mesh_->cell_map(false),1);
  
  int num_with_cells = 0;
  for (const_op_iterator it = OpBegin(); it != OpEnd(); ++it) {
    if ((*it)->schema & OPERATOR_SCHEMA_DOFS_CELL) {
      if ((*it)->schema == (OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_CELL)) {
        // diagonal schema
        for (int c = 0; c != ncells_owned; ++c) {
          D_c[0][c] += (*it)->vals[c];
        }
      } else {
        num_with_cells++;
      }
    }
  }
  if (num_with_cells == 0 || num_with_cells > 1) {
    Errors::Message msg("Schur complement to Sff must have exactly one cell-based schema with off-diagonal entries or other, non-cell dofs.");
    Exceptions::amanzi_throw(msg);
  }
  
  // schur complement
  int i_schur = 0;
  for (const_op_iterator it = OpBegin(); it != OpEnd(); ++it) {
    if ((*it)->schema == (OPERATOR_SCHEMA_BASE_CELL |
                          OPERATOR_SCHEMA_DOFS_CELL | OPERATOR_SCHEMA_DOFS_FACE)) {
      ASSERT((*it)->matrices.size() == ncells_owned);

      // create or get extra ops, and keep them for future use
      const std::vector<int>& bc_model = bc_->bc_model();
      Teuchos::RCP<Op_Cell_Face> schur_op;

      if (schur_ops_.size() > i_schur) {
        // get existing
        schur_op = schur_ops_[i_schur];
      } else {
        // create and fill
        ASSERT(i_schur == schur_ops_.size());
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
        //        if (tmp == 0.0 && (*K_)[c].isZero()) continue;  // We skip zero matrices

        for (int n = 0; n < nfaces; n++) {
          for (int m = 0; m < nfaces; m++) {
            Scell(n, m) = Acell(n, m) - Acell(n, nfaces) * Acell(nfaces, m) / tmp;
          }
        }

        for (int n = 0; n < nfaces; n++) {  // Symbolic boundary conditions
          int f = faces[n];
          if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
            for (int m = 0; m < nfaces; m++) Scell(n, m) = Scell(m, n) = 0.0;
            Scell(n, n) = 1.0;
          }
        }
        ASSERT(true);
      }

      // Assemble this Schur Op into matrix
      schur_op->AssembleMatrixOp(this, map, matrix, my_block_row, my_block_col);
    }
  }
}

// visit method for Apply -- this is identical to Operator_FaceCell's
// version.
int
Operator_FaceCellSff::ApplyMatrixFreeOp(const Op_Cell_FaceCell& op,
        const CompositeVector& X, CompositeVector& Y) const
{
  ASSERT(op.matrices.size() == ncells_owned);

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
void
Operator_FaceCellSff::SymbolicAssembleMatrix()
{
  // SuperMap for Sff is face only
  CompositeVectorSpace smap_space;
  smap_space.SetMesh(cvs_->Mesh())->SetComponent("face", AmanziMesh::FACE, 1);
  smap_ = createSuperMap(smap_space, schema(), 1);

  // create the graph
  int row_size = MaxRowSize(*mesh_, schema(), 1);
  Teuchos::RCP<GraphFE> graph = Teuchos::rcp(new GraphFE(smap_->Map(),
          smap_->GhostedMap(), smap_->GhostedMap(), row_size));

  Operator::SymbolicAssembleMatrix(*smap_, *graph, 0, 0);

  // Completing and optimizing the graphs
  int ierr = graph->FillComplete(smap_->Map(), smap_->Map());
  ASSERT(!ierr);

  // create global matrix
  Amat_ = Teuchos::rcp(new MatrixFE(graph));
  A_ = Amat_->Matrix();
}


// visit method for sparsity structure of Schur complement
void
Operator_FaceCellSff::SymbolicAssembleMatrixOp(const Op_Cell_FaceCell& op,
        const SuperMap& map, GraphFE& graph,
        int my_block_row, int my_block_col) const
{
  std::string name = "Sff alt CELL_FACE";
  Op_Cell_Face schur_op(name, mesh_);
  Operator_FaceCell::SymbolicAssembleMatrixOp(schur_op, map, graph, my_block_row, my_block_col);
}


}  // namespace Operators
}  // namespace Amanzi

