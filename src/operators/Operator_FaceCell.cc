/*
  This is the Operator component of the Amanzi code.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
  Ethan Coon (ecoon@lanl.gov)
*/

#include "DenseMatrix.hh"
#include "Op_Cell_FaceCell.hh"
#include "Op_Cell_Face.hh"

#include "SuperMap.hh"
#include "GraphFE.hh"
#include "MatrixFE.hh"

#include "OperatorDefs.hh"
#include "Operator_FaceCell.hh"

/* ******************************************************************
   Operator whose unknowns are CELL + FACE

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

// visit methods for Apply

/* ******************************************************************
* Apply the local matrices directly as schemas match.
****************************************************************** */
int
Operator_FaceCell::ApplyMatrixFreeOp(const Op_Cell_FaceCell& op,
        const CompositeVector& X, CompositeVector& Y) const
{
  ASSERT(op.matrices.size() == ncells_owned);

  X.ScatterMasterToGhosted();
  const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell");

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


int
Operator_FaceCell::ApplyMatrixFreeOp(const Op_Cell_Face& op,
        const CompositeVector& X, CompositeVector& Y) const
{
  ASSERT(op.matrices.size() == ncells_owned);

  X.ScatterMasterToGhosted();
  const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

  {
    Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);

    AmanziMesh::Entity_ID_List faces;
    for (int c=0; c!=ncells_owned; ++c) {
      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      WhetStone::DenseVector v(nfaces), av(nfaces);
      for (int n=0; n!=nfaces; ++n) {
        v(n) = Xf[0][faces[n]];
      }

      const WhetStone::DenseMatrix& Acell = op.matrices[c];
      Acell.Multiply(v, av, false);

      for (int n=0; n!=nfaces; ++n) {
        Yf[0][faces[n]] += av(n);
      }
    } 
  }
  Y.GatherGhostedToMaster(Add);
  return 0;
}


// visit methods for symbolic assemble
void
Operator_FaceCell::SymbolicAssembleMatrixOp(const Op_Cell_FaceCell& op,
        const SuperMap& map, GraphFE& graph,
        int my_block_row, int my_block_col) const
{
  int lid_r[OPERATOR_MAX_FACES];
  int lid_c[OPERATOR_MAX_FACES];

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices("face", my_block_row);
  const std::vector<int>& face_col_inds = map.GhostIndices("face", my_block_col);
  const std::vector<int>& cell_row_inds = map.GhostIndices("cell", my_block_row);
  const std::vector<int>& cell_col_inds = map.GhostIndices("cell", my_block_col);

  int ierr(0);
  AmanziMesh::Entity_ID_List faces;
  for (int c=0; c!=ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int n=0; n!=nfaces; ++n) {
      lid_r[n] = face_row_inds[faces[n]];
      lid_c[n] = face_col_inds[faces[n]];
    }
    lid_r[nfaces] = cell_row_inds[c];
    lid_c[nfaces] = cell_col_inds[c];
    ierr |= graph.InsertMyIndices(nfaces+1, lid_r, nfaces+1, lid_c);
  }
  ASSERT(!ierr);
}

void
Operator_FaceCell::SymbolicAssembleMatrixOp(const Op_Cell_Face& op,
        const SuperMap& map, GraphFE& graph,
        int my_block_row, int my_block_col) const
{
  int lid_r[OPERATOR_MAX_FACES];
  int lid_c[OPERATOR_MAX_FACES];

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices("face", my_block_row);
  const std::vector<int>& face_col_inds = map.GhostIndices("face", my_block_col);

  int ierr(0);
  AmanziMesh::Entity_ID_List faces;
  for (int c=0; c!=ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int n=0; n!=nfaces; ++n) {
      lid_r[n] = face_row_inds[faces[n]];
      lid_c[n] = face_col_inds[faces[n]];
    }
    ierr |= graph.InsertMyIndices(nfaces, lid_r, nfaces, lid_c);
  }
  ASSERT(!ierr);
}

// visit methods for assemble
void
Operator_FaceCell::AssembleMatrixOp(const Op_Cell_FaceCell& op,
        const SuperMap& map, MatrixFE& mat,
        int my_block_row, int my_block_col) const
{
  ASSERT(op.matrices.size() == ncells_owned);

  int lid_r[OPERATOR_MAX_FACES+1];
  int lid_c[OPERATOR_MAX_FACES+1];

  // ELEMENT: cell, DOFS: face and cell
  const std::vector<int>& face_row_inds = map.GhostIndices("face", my_block_row);
  const std::vector<int>& face_col_inds = map.GhostIndices("face", my_block_col);
  const std::vector<int>& cell_row_inds = map.GhostIndices("cell", my_block_row);
  const std::vector<int>& cell_col_inds = map.GhostIndices("cell", my_block_col);        

  int ierr(0);
  AmanziMesh::Entity_ID_List faces;
  for (int c=0; c!=ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    
    int nfaces = faces.size();
    for (int n=0; n!=nfaces; ++n) {
      lid_r[n] = face_row_inds[faces[n]];
      lid_c[n] = face_col_inds[faces[n]];
    }
    lid_r[nfaces] = cell_row_inds[c];
    lid_c[nfaces] = cell_col_inds[c];

    ierr |= mat.SumIntoMyValues(lid_r, lid_c, op.matrices[c]);
  }
  ASSERT(!ierr);
}


void
Operator_FaceCell::AssembleMatrixOp(const Op_Cell_Face& op,
        const SuperMap& map, MatrixFE& mat,
        int my_block_row, int my_block_col) const
{
  ASSERT(op.matrices.size() == ncells_owned);

  int lid_r[OPERATOR_MAX_FACES];
  int lid_c[OPERATOR_MAX_FACES];

  // ELEMENT: cell, DOFS: face and cell
  const std::vector<int>& face_row_inds = map.GhostIndices("face", my_block_row);
  const std::vector<int>& face_col_inds = map.GhostIndices("face", my_block_col);

  int ierr(0);
  AmanziMesh::Entity_ID_List faces;
  for (int c=0; c!=ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    
    int nfaces = faces.size();
    for (int n=0; n!=nfaces; ++n) {
      lid_r[n] = face_row_inds[faces[n]];
      lid_c[n] = face_col_inds[faces[n]];
    }
    
    ierr |= mat.SumIntoMyValues(lid_r, lid_c, op.matrices[c]);
  }
  ASSERT(!ierr);
}


}  // namespace Operators
}  // namespace Amanzi



