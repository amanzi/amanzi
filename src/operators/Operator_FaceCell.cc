/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  Operator whose unknowns are CELL + FACE.
*/

#include "DenseMatrix.hh"
#include "Op_Cell_FaceCell.hh"
#include "Op_Cell_Face.hh"
#include "Op_SurfaceCell_SurfaceCell.hh"
#include "Op_SurfaceFace_SurfaceCell.hh"

#include "SuperMap.hh"
#include "GraphFE.hh"
#include "MatrixFE.hh"

#include "OperatorDefs.hh"
#include "Operator_FaceCell.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Visit methods for Apply:
* apply the local matrices directly as schemas match.
****************************************************************** */
int Operator_FaceCell::ApplyMatrixFreeOp(const Op_Cell_FaceCell& op,
                                         const CompositeVector& X, CompositeVector& Y) const
{
  ASSERT(op.matrices.size() == ncells_owned);

  Y.PutScalarGhosted(0.);
  X.ScatterMasterToGhosted();
  const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell");

  {
    Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);
    Epetra_MultiVector& Yc = *Y.ViewComponent("cell");

    AmanziMesh::Entity_ID_List faces;
    for (int c = 0; c != ncells_owned; ++c) {
      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      WhetStone::DenseVector v(nfaces + 1), av(nfaces + 1);
      for (int n = 0; n != nfaces; ++n) {
        v(n) = Xf[0][faces[n]];
      }
      v(nfaces) = Xc[0][c];

      const WhetStone::DenseMatrix& Acell = op.matrices[c];
      Acell.Multiply(v, av, false);

      for (int n = 0; n != nfaces; ++n) {
        Yf[0][faces[n]] += av(n);
      }
      Yc[0][c] += av(nfaces);
    } 
  }

  Y.GatherGhostedToMaster(Add);
  return 0;
}


/* ******************************************************************
* Apply the local matrices directly as schemas match.
****************************************************************** */
int Operator_FaceCell::ApplyMatrixFreeOp(const Op_Cell_Face& op,
                                         const CompositeVector& X, CompositeVector& Y) const
{
  ASSERT(op.matrices.size() == ncells_owned);

  Y.PutScalarGhosted(0.);
  X.ScatterMasterToGhosted();
  const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

  {
    Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);

    AmanziMesh::Entity_ID_List faces;
    for (int c = 0; c != ncells_owned; ++c) {
      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      WhetStone::DenseVector v(nfaces), av(nfaces);
      for (int n = 0; n != nfaces; ++n) {
        v(n) = Xf[0][faces[n]];
      }

      const WhetStone::DenseMatrix& Acell = op.matrices[c];
      Acell.Multiply(v, av, false);

      for (int n = 0; n != nfaces; ++n) {
        Yf[0][faces[n]] += av(n);
      }
    } 
  }

  Y.GatherGhostedToMaster(Add);
  return 0;
}


/* ******************************************************************
* visit method for apply surface cells into subsurface faces
****************************************************************** */
int Operator_FaceCell::ApplyMatrixFreeOp(const Op_SurfaceCell_SurfaceCell& op,
                                         const CompositeVector& X, CompositeVector& Y) const
{
  int nsurf_cells = op.surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  ASSERT(op.vals.size() == nsurf_cells);

  const Epetra_MultiVector& Xf = *X.ViewComponent("face", false);
  Epetra_MultiVector& Yf = *Y.ViewComponent("face", false);
  for (int sc = 0; sc != nsurf_cells; ++sc) {
    int f = op.surf_mesh->entity_get_parent(AmanziMesh::CELL, sc);
    Yf[0][f] += op.vals[sc] * Xf[0][f];
  } 
  return 0;
}


/* ******************************************************************
* visit method for apply surface cells into subsurface faces
****************************************************************** */
int Operator_FaceCell::ApplyMatrixFreeOp(const Op_SurfaceFace_SurfaceCell& op,
                                         const CompositeVector& X, CompositeVector& Y) const
{
  int nsurf_faces = op.surf_mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  ASSERT(op.matrices.size() == nsurf_faces);

  X.ScatterMasterToGhosted();
  const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

  Y.PutScalarGhosted(0.);
  Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);
  
  AmanziMesh::Entity_ID_List cells;
  for (int sf = 0; sf != nsurf_faces; ++sf) {
    op.surf_mesh->face_get_cells(sf, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    WhetStone::DenseVector v(ncells), av(ncells);
    for (int n = 0; n != ncells; ++n) {
      v(n) = Xf[0][op.surf_mesh->entity_get_parent(AmanziMesh::CELL, cells[n])];
    }

    const WhetStone::DenseMatrix& Aface = op.matrices[sf];
    Aface.Multiply(v, av, false);

    for (int n = 0; n != ncells; ++n) {
      Yf[0][op.surf_mesh->entity_get_parent(AmanziMesh::CELL, cells[n])] += av(n);
    }
  } 
  Y.GatherGhostedToMaster("face");
  return 0;
}


/* ******************************************************************
* Visit methods for symbolic assemble: FaceCell
****************************************************************** */
void Operator_FaceCell::SymbolicAssembleMatrixOp(const Op_Cell_FaceCell& op,
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
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int n = 0; n != nfaces; ++n) {
      lid_r[n] = face_row_inds[faces[n]];
      lid_c[n] = face_col_inds[faces[n]];
    }
    lid_r[nfaces] = cell_row_inds[c];
    lid_c[nfaces] = cell_col_inds[c];
    ierr |= graph.InsertMyIndices(nfaces+1, lid_r, nfaces+1, lid_c);
  }
  ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for symbolic assemble: Face
****************************************************************** */
void Operator_FaceCell::SymbolicAssembleMatrixOp(const Op_Cell_Face& op,
                                                 const SuperMap& map, GraphFE& graph,
                                                 int my_block_row, int my_block_col) const
{
  int lid_r[OPERATOR_MAX_FACES];
  int lid_c[OPERATOR_MAX_FACES];

  // ELEMENT: cell, DOFS: face
  const std::vector<int>& face_row_inds = map.GhostIndices("face", my_block_row);
  const std::vector<int>& face_col_inds = map.GhostIndices("face", my_block_col);

  int ierr(0);
  AmanziMesh::Entity_ID_List faces;
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int n = 0; n != nfaces; ++n) {
      lid_r[n] = face_row_inds[faces[n]];
      lid_c[n] = face_col_inds[faces[n]];
    }
    ierr |= graph.InsertMyIndices(nfaces, lid_r, nfaces, lid_c);
  }
  ASSERT(!ierr);
}


void
Operator_FaceCell::SymbolicAssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op,
        const SuperMap& map, GraphFE& graph,
        int my_block_row, int my_block_col) const
{
  int nsurf_cells = op.surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices("face", my_block_row);
  const std::vector<int>& face_col_inds = map.GhostIndices("face", my_block_col);

  int ierr = 0;
  for (int sc=0; sc!=nsurf_cells; ++sc) {
    int lid_r = face_row_inds[op.surf_mesh->entity_get_parent(AmanziMesh::CELL,sc)];
    int lid_c = face_col_inds[op.surf_mesh->entity_get_parent(AmanziMesh::CELL,sc)];
    ierr |= graph.InsertMyIndices(lid_r, 1, &lid_c);
  }
  ASSERT(!ierr);
}


void
Operator_FaceCell::SymbolicAssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op,
        const SuperMap& map, GraphFE& graph,
        int my_block_row, int my_block_col) const
{
  int nsurf_faces = op.surf_mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int lid_r[2];
  int lid_c[2];

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices("face", my_block_row);
  const std::vector<int>& face_col_inds = map.GhostIndices("face", my_block_col);

  int ierr = 0;
  AmanziMesh::Entity_ID_List cells;
  for (int sf=0; sf!=nsurf_faces; ++sf) {
    op.surf_mesh->face_get_cells(sf, AmanziMesh::USED, &cells);
    int ncells = cells.size();
    for (int n=0; n!=ncells; ++n) {
      int f = op.surf_mesh->entity_get_parent(AmanziMesh::CELL,cells[n]);
      lid_r[n] = face_row_inds[f];
      lid_c[n] = face_col_inds[f];
    }

    ierr |= graph.InsertMyIndices(ncells, lid_r, ncells, lid_c);
  }
  ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble: FaceCell
****************************************************************** */
void Operator_FaceCell::AssembleMatrixOp(const Op_Cell_FaceCell& op,
                                         const SuperMap& map, MatrixFE& mat,
                                         int my_block_row, int my_block_col) const
{
  ASSERT(op.matrices.size() == ncells_owned);

  int lid_r[OPERATOR_MAX_FACES + 1];
  int lid_c[OPERATOR_MAX_FACES + 1];

  // ELEMENT: cell, DOFS: face and cell
  const std::vector<int>& face_row_inds = map.GhostIndices("face", my_block_row);
  const std::vector<int>& face_col_inds = map.GhostIndices("face", my_block_col);
  const std::vector<int>& cell_row_inds = map.GhostIndices("cell", my_block_row);
  const std::vector<int>& cell_col_inds = map.GhostIndices("cell", my_block_col);        

  int ierr(0);
  AmanziMesh::Entity_ID_List faces;
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    
    int nfaces = faces.size();
    for (int n = 0; n != nfaces; ++n) {
      lid_r[n] = face_row_inds[faces[n]];
      lid_c[n] = face_col_inds[faces[n]];
    }
    lid_r[nfaces] = cell_row_inds[c];
    lid_c[nfaces] = cell_col_inds[c];

    ierr |= mat.SumIntoMyValues(lid_r, lid_c, op.matrices[c]);
  }
  ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble: Face
****************************************************************** */
void Operator_FaceCell::AssembleMatrixOp(const Op_Cell_Face& op,
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
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    
    int nfaces = faces.size();
    for (int n = 0; n != nfaces; ++n) {
      lid_r[n] = face_row_inds[faces[n]];
      lid_c[n] = face_col_inds[faces[n]];
    }
    
    ierr |= mat.SumIntoMyValues(lid_r, lid_c, op.matrices[c]);
  }
  ASSERT(!ierr);
}


void
Operator_FaceCell::AssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op,
        const SuperMap& map, MatrixFE& mat,
        int my_block_row, int my_block_col) const
{
  int nsurf_cells = op.surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices("face", my_block_row);
  const std::vector<int>& face_col_inds = map.GhostIndices("face", my_block_col);

  int ierr = 0;
  for (int sc=0; sc!=nsurf_cells; ++sc) {
    int lid_r = face_row_inds[op.surf_mesh->entity_get_parent(AmanziMesh::CELL,sc)];
    int lid_c = face_col_inds[op.surf_mesh->entity_get_parent(AmanziMesh::CELL,sc)];
    ierr |= mat.SumIntoMyValues(lid_r, 1, &op.vals[sc], &lid_c);
  }
  ASSERT(!ierr);
}


void
Operator_FaceCell::AssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op,
        const SuperMap& map, MatrixFE& mat,
        int my_block_row, int my_block_col) const
{
  int nsurf_faces = op.surf_mesh->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  int lid_r[2];
  int lid_c[2];

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices("face", my_block_row);
  const std::vector<int>& face_col_inds = map.GhostIndices("face", my_block_col);

  int ierr = 0;
  AmanziMesh::Entity_ID_List cells;
  for (int sf=0; sf!=nsurf_faces; ++sf) {
    op.surf_mesh->face_get_cells(sf, AmanziMesh::USED, &cells);
    int ncells = cells.size();
    for (int n=0; n!=ncells; ++n) {
      int f = op.surf_mesh->entity_get_parent(AmanziMesh::CELL,cells[n]);
      lid_r[n] = face_row_inds[f];
      lid_c[n] = face_col_inds[f];
    }

    ierr |= mat.SumIntoMyValues(lid_r, lid_c, op.matrices[sf]);
  }
  ASSERT(!ierr);
}

}  // namespace Operators
}  // namespace Amanzi



