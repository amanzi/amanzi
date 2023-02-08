/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)
*/

/*
  Operators

  Operator whose unknowns are CELL + FACE.
*/

#include "DenseMatrix.hh"
#include "Op_Cell_FaceCell.hh"
#include "Op_Cell_Face.hh"
#include "Op_Diagonal.hh"
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
int
Operator_FaceCell::ApplyMatrixFreeOp(const Op_Cell_FaceCell& op,
                                     const CompositeVector& X,
                                     CompositeVector& Y) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);
  const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell");

  {
    Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);
    Epetra_MultiVector& Yc = *Y.ViewComponent("cell");

    const auto& map = Yf.Map();

    for (int c = 0; c != ncells_owned; ++c) {
      const auto& faces = mesh_->getCellFaces(c);
      int nfaces = faces.size();

      int npoints(0);
      for (int n = 0; n != nfaces; ++n) npoints += map.ElementSize(faces[n]);

      int m(0);
      WhetStone::DenseVector v(npoints + 1), av(npoints + 1);

      for (int n = 0; n != nfaces; ++n) {
        int f = faces[n];
        int first = map.FirstPointInElement(f);
        for (int k = 0; k < map.ElementSize(f); ++k) { v(m++) = Xf[0][first + k]; }
      }
      v(npoints) = Xc[0][c];

      const WhetStone::DenseMatrix& Acell = op.matrices[c];
      Acell.Multiply(v, av, false);

      m = 0;
      for (int n = 0; n != nfaces; ++n) {
        int f = faces[n];
        int first = map.FirstPointInElement(f);
        for (int k = 0; k < map.ElementSize(f); ++k) { Yf[0][first + k] += av(m++); }
      }
      Yc[0][c] += av(npoints);
    }
  }
  return 0;
}


/* ******************************************************************
* Apply the local matrices directly as schemas match.
****************************************************************** */
int
Operator_FaceCell::ApplyMatrixFreeOp(const Op_Cell_Face& op,
                                     const CompositeVector& X,
                                     CompositeVector& Y) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);
  const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

  {
    Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);

    for (int c = 0; c != ncells_owned; ++c) {
      const auto& faces = mesh_->getCellFaces(c);
      int nfaces = faces.size();

      WhetStone::DenseVector v(nfaces), av(nfaces);
      for (int n = 0; n != nfaces; ++n) { v(n) = Xf[0][faces[n]]; }

      const WhetStone::DenseMatrix& Acell = op.matrices[c];
      Acell.Multiply(v, av, false);

      for (int n = 0; n != nfaces; ++n) { Yf[0][faces[n]] += av(n); }
    }
  }
  return 0;
}


/* ******************************************************************
* visit method for apply surface cells into subsurface faces
****************************************************************** */
int
Operator_FaceCell::ApplyMatrixFreeOp(const Op_SurfaceCell_SurfaceCell& op,
                                     const CompositeVector& X,
                                     CompositeVector& Y) const
{
  int nsurf_cells = op.surf_mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  AMANZI_ASSERT(op.diag->MyLength() == nsurf_cells);

  const Epetra_MultiVector& Xf = *X.ViewComponent("face", false);
  Epetra_MultiVector& Yf = *Y.ViewComponent("face", false);
  for (int sc = 0; sc != nsurf_cells; ++sc) {
    int f = op.surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
    Yf[0][f] += (*op.diag)[0][sc] * Xf[0][f];
  }
  return 0;
}


/* ******************************************************************
* visit method for apply surface cells into subsurface faces
****************************************************************** */
int
Operator_FaceCell::ApplyMatrixFreeOp(const Op_SurfaceFace_SurfaceCell& op,
                                     const CompositeVector& X,
                                     CompositeVector& Y) const
{
  int nsurf_faces = op.surf_mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  AMANZI_ASSERT(op.matrices.size() == nsurf_faces);
  const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);
  Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);

  AmanziMesh::Entity_ID_View cells;
  for (int sf = 0; sf != nsurf_faces; ++sf) {
    cells = op.surf_mesh->getFaceCells(sf, AmanziMesh::Parallel_kind::ALL);
    int ncells = cells.size();

    WhetStone::DenseVector v(ncells), av(ncells);
    for (int n = 0; n != ncells; ++n) {
      v(n) = Xf[0][op.surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, cells[n])];
    }

    const WhetStone::DenseMatrix& Aface = op.matrices[sf];
    Aface.Multiply(v, av, false);

    for (int n = 0; n != ncells; ++n) {
      Yf[0][op.surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, cells[n])] += av(n);
    }
  }
  return 0;
}


/* ******************************************************************
* Visit methods for symbolic assemble: FaceCell
****************************************************************** */
void
Operator_FaceCell::SymbolicAssembleMatrixOp(const Op_Cell_FaceCell& op,
                                            const SuperMap& map,
                                            GraphFE& graph,
                                            int my_block_row,
                                            int my_block_col) const
{
  std::vector<int> lid_r(2 * cell_max_faces + 1);
  std::vector<int> lid_c(2 * cell_max_faces + 1);

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const std::vector<int>& face_col_inds = map.GhostIndices(my_block_col, "face", 0);
  const std::vector<int>& cell_row_inds = map.GhostIndices(my_block_row, "cell", 0);
  const std::vector<int>& cell_col_inds = map.GhostIndices(my_block_col, "cell", 0);

  Teuchos::RCP<const Epetra_BlockMap> face_gh_map = map.ComponentGhostedMap(my_block_row, "face");
  Teuchos::RCP<const Epetra_BlockMap> cell_map = map.ComponentMap(my_block_row, "cell");

  int ierr(0);
  for (int c = 0; c != ncells_owned; ++c) {
    const auto& faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    int k = 0;
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      int face_dof_size = face_gh_map->ElementSize(f);
      int first = face_gh_map->FirstPointInElement(f);

      for (int m = 0; m != face_dof_size; ++m) {
        lid_r[k] = face_row_inds[first + m];
        lid_c[k] = face_col_inds[first + m];
        k++;
      }
    }

    int cell_dof_size = cell_map->ElementSize(c);
    int first = cell_map->FirstPointInElement(c);
    for (int m = 0; m != cell_dof_size; ++m) {
      lid_r[k] = cell_row_inds[first + m];
      lid_c[k] = cell_col_inds[first + m];
      k++;
    }

    ierr |= graph.InsertMyIndices(k, lid_r.data(), k, lid_c.data());
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for symbolic assemble: Face
****************************************************************** */
void
Operator_FaceCell::SymbolicAssembleMatrixOp(const Op_Cell_Face& op,
                                            const SuperMap& map,
                                            GraphFE& graph,
                                            int my_block_row,
                                            int my_block_col) const
{
  std::vector<int> lid_r(2 * cell_max_faces);
  std::vector<int> lid_c(2 * cell_max_faces);

  // ELEMENT: cell, DOFS: face
  const std::vector<int>& face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const std::vector<int>& face_col_inds = map.GhostIndices(my_block_col, "face", 0);

  Teuchos::RCP<const Epetra_BlockMap> face_gh_map = map.ComponentGhostedMap(my_block_row, "face");

  int ierr(0);
  for (int c = 0; c != ncells_owned; ++c) {
    const auto& faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    int k = 0;
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      int face_dof_size = face_gh_map->ElementSize(f);
      int first = face_gh_map->FirstPointInElement(f);

      for (int m = 0; m != face_dof_size; ++m) {
        lid_r[k] = face_row_inds[first + m];
        lid_c[k] = face_col_inds[first + m];
        k++;
      }
    }
    ierr |= graph.InsertMyIndices(k, lid_r.data(), k, lid_c.data());
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for symbolic assemble: Surface
****************************************************************** */
void
Operator_FaceCell::SymbolicAssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op,
                                            const SuperMap& map,
                                            GraphFE& graph,
                                            int my_block_row,
                                            int my_block_col) const
{
  int nsurf_cells = op.surf_mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const std::vector<int>& face_col_inds = map.GhostIndices(my_block_col, "face", 0);

  int ierr = 0;
  for (int sc = 0; sc != nsurf_cells; ++sc) {
    int lid_r = face_row_inds[op.surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, sc)];
    int lid_c = face_col_inds[op.surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, sc)];
    ierr |= graph.InsertMyIndices(lid_r, 1, &lid_c);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for symbolic assemble: Surface
****************************************************************** */
void
Operator_FaceCell::SymbolicAssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op,
                                            const SuperMap& map,
                                            GraphFE& graph,
                                            int my_block_row,
                                            int my_block_col) const
{
  int nsurf_faces = op.surf_mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  int lid_r[2];
  int lid_c[2];

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const std::vector<int>& face_col_inds = map.GhostIndices(my_block_col, "face", 0);

  int ierr = 0;
  AmanziMesh::Entity_ID_View cells;
  for (int sf = 0; sf != nsurf_faces; ++sf) {
    cells = op.surf_mesh->getFaceCells(sf, AmanziMesh::Parallel_kind::ALL);
    int ncells = cells.size();
    for (int n = 0; n != ncells; ++n) {
      int f = op.surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, cells[n]);
      lid_r[n] = face_row_inds[f];
      lid_c[n] = face_col_inds[f];
    }

    ierr |= graph.InsertMyIndices(ncells, lid_r, ncells, lid_c);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble: FaceCell
****************************************************************** */
void
Operator_FaceCell::AssembleMatrixOp(const Op_Cell_FaceCell& op,
                                    const SuperMap& map,
                                    MatrixFE& mat,
                                    int my_block_row,
                                    int my_block_col) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

  std::vector<int> lid_r(2 * cell_max_faces + 1);
  std::vector<int> lid_c(2 * cell_max_faces + 1);

  // ELEMENT: cell, DOFS: face and cell
  const std::vector<int>& face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const std::vector<int>& face_col_inds = map.GhostIndices(my_block_col, "face", 0);
  const std::vector<int>& cell_row_inds = map.GhostIndices(my_block_row, "cell", 0);
  const std::vector<int>& cell_col_inds = map.GhostIndices(my_block_col, "cell", 0);

  Teuchos::RCP<const Epetra_BlockMap> face_gh_map = map.ComponentGhostedMap(my_block_row, "face");
  Teuchos::RCP<const Epetra_BlockMap> cell_map = map.ComponentMap(my_block_row, "cell");

  int ierr(0);
  for (int c = 0; c != ncells_owned; ++c) {
    const auto& faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    int k = 0;
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      int face_dof_size = face_gh_map->ElementSize(f);
      int first = face_gh_map->FirstPointInElement(f);

      for (int m = 0; m != face_dof_size; ++m) {
        lid_r[k] = face_row_inds[first + m];
        lid_c[k] = face_col_inds[first + m];
        k++;
      }
    }

    int cell_dof_size = cell_map->ElementSize(c);
    int first = cell_map->FirstPointInElement(c);
    for (int m = 0; m != cell_dof_size; ++m) {
      lid_r[k] = cell_row_inds[first + m];
      lid_c[k] = cell_col_inds[first + m];
      k++;
    }

    ierr |= mat.SumIntoMyValues(lid_r.data(), lid_c.data(), op.matrices[c]);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble: Face
****************************************************************** */
void
Operator_FaceCell::AssembleMatrixOp(const Op_Cell_Face& op,
                                    const SuperMap& map,
                                    MatrixFE& mat,
                                    int my_block_row,
                                    int my_block_col) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

  std::vector<int> lid_r(cell_max_faces + 1);
  std::vector<int> lid_c(cell_max_faces + 1);

  // ELEMENT: cell, DOFS: face and cell
  const std::vector<int>& face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const std::vector<int>& face_col_inds = map.GhostIndices(my_block_col, "face", 0);

  Teuchos::RCP<const Epetra_BlockMap> face_gh_map = map.ComponentGhostedMap(my_block_row, "face");

  int ierr(0);
  for (int c = 0; c != ncells_owned; ++c) {
    const auto& faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    int k = 0;
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      int face_dof_size = face_gh_map->ElementSize(f);
      int first = face_gh_map->FirstPointInElement(f);

      for (int m = 0; m != face_dof_size; ++m) {
        lid_r[k] = face_row_inds[first + m];
        lid_c[k] = face_col_inds[first + m];
        k++;
      }
    }

    ierr |= mat.SumIntoMyValues(lid_r.data(), lid_c.data(), op.matrices[c]);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble: Surface
****************************************************************** */
void
Operator_FaceCell::AssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op,
                                    const SuperMap& map,
                                    MatrixFE& mat,
                                    int my_block_row,
                                    int my_block_col) const
{
  int nsurf_cells = op.surf_mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const std::vector<int>& face_col_inds = map.GhostIndices(my_block_col, "face", 0);

  int ierr = 0;
  for (int sc = 0; sc != nsurf_cells; ++sc) {
    int lid_r = face_row_inds[op.surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, sc)];
    int lid_c = face_col_inds[op.surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, sc)];
    ierr |= mat.SumIntoMyValues(lid_r, 1, &(*op.diag)[0][sc], &lid_c);
  }
  AMANZI_ASSERT(!ierr);
}


void
Operator_FaceCell::AssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op,
                                    const SuperMap& map,
                                    MatrixFE& mat,
                                    int my_block_row,
                                    int my_block_col) const
{
  int nsurf_faces = op.surf_mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  int lid_r[2];
  int lid_c[2];

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const std::vector<int>& face_col_inds = map.GhostIndices(my_block_col, "face", 0);

  int ierr = 0;
  AmanziMesh::Entity_ID_View cells;
  for (int sf = 0; sf != nsurf_faces; ++sf) {
    cells = op.surf_mesh->getFaceCells(sf, AmanziMesh::Parallel_kind::ALL);
    int ncells = cells.size();
    for (int n = 0; n != ncells; ++n) {
      int f = op.surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, cells[n]);
      lid_r[n] = face_row_inds[f];
      lid_c[n] = face_col_inds[f];
    }

    ierr |= mat.SumIntoMyValues(lid_r, lid_c, op.matrices[sf]);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Operator>
Operator_FaceCell::Clone() const
{
  return Teuchos::rcp(new Operator_FaceCell(*this));
}

} // namespace Operators
} // namespace Amanzi
