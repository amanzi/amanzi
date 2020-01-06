/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  Operator whose unknowns are CELL
*/

#include "DenseMatrix.hh"
#include "Op_Cell_FaceCell.hh"

#include "SuperMap.hh"
#include "GraphFE.hh"
#include "MatrixFE.hh"
#include "Operator_ConsistentFace.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Visit methods for Apply.
* Apply the local matrices directly as schema is a subset of 
* assembled schema.
****************************************************************** */
int Operator_ConsistentFace::ApplyMatrixFreeOp(
    const Op_Cell_FaceCell& op, const CompositeVector& X, CompositeVector& Y) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

  Y.PutScalarGhosted(0.);
  X.ScatterMasterToGhosted();
  const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

  {
    Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);

    AmanziMesh::Entity_ID_List faces;
    for (int c=0; c!=ncells_owned; ++c) {
      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      WhetStone::DenseVector v(nfaces), av(nfaces);
      av.PutScalar(0.0);
      for (int n=0; n!=nfaces; ++n) {
        v(n) = Xf[0][faces[n]];
      }

      const WhetStone::DenseMatrix& Acell = op.matrices[c];
      // must do multiply manually because Acell is not nfaces x nfaces
      for (int n=0; n!=nfaces; ++n)
        for (int m=0; m!=nfaces; ++m)
          av(m) += Acell(m,n) * v(n);

      for (int n=0; n!=nfaces; ++n) {
        Yf[0][faces[n]] += av(n);
        AMANZI_ASSERT(std::abs(av(n)) < 1.e20);
      }
    } 
  }
  Y.GatherGhostedToMaster("face", Add);
  return 0;
}


/* ******************************************************************
* Visit methods for symbolic assemble.
* Insert the diagonal on cells
****************************************************************** */
void Operator_ConsistentFace::SymbolicAssembleMatrixOp(
    const Op_Cell_FaceCell& op,
    const SuperMap& map, GraphFE& graph,
    int my_block_row, int my_block_col) const
{
  std::vector<int> lid_r(cell_max_faces);
  std::vector<int> lid_c(cell_max_faces);

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const std::vector<int>& face_col_inds = map.GhostIndices(my_block_col, "face", 0);

  int ierr(0);
  AmanziMesh::Entity_ID_List faces;
  for (int c=0; c!=ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int n=0; n!=nfaces; ++n) {
      lid_r[n] = face_row_inds[faces[n]];
      lid_c[n] = face_col_inds[faces[n]];
    }
    ierr |= graph.InsertMyIndices(nfaces, lid_r.data(), nfaces, lid_c.data());
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble
* Insert each cells neighboring cells.
****************************************************************** */
void Operator_ConsistentFace::AssembleMatrixOp(
    const Op_Cell_FaceCell& op,
    const SuperMap& map, MatrixFE& mat,
    int my_block_row, int my_block_col) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

  std::vector<int> lid_r(cell_max_faces);
  std::vector<int> lid_c(cell_max_faces);
  std::vector<double> vals(cell_max_faces);

  // ELEMENT: cell, DOFS: face and cell
  const std::vector<int>& face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const std::vector<int>& face_col_inds = map.GhostIndices(my_block_col, "face", 0);

  int ierr(0);
  AmanziMesh::Entity_ID_List faces;
  for (int c=0; c!=ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    
    int nfaces = faces.size();
    for (int n=0; n!=nfaces; ++n) {
      lid_r[n] = face_row_inds[faces[n]];
      lid_c[n] = face_col_inds[faces[n]];
    }

    for (int n=0; n!=nfaces; ++n) {
      for (int m=0; m!=nfaces; ++m) vals[m] = op.matrices[c](n,m);
      ierr |= mat.SumIntoMyValues(lid_r[n], nfaces, vals.data(), lid_c.data());
    }
  }
  AMANZI_ASSERT(!ierr);
}

}  // namespace Operators
}  // namespace Amanzi

