/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

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
  AMANZI_ASSERT(op.data.extent(0) == ncells_owned);

  Y.putScalarGhosted(0.0);
  X.ScatterMasterToGhosted();
  const auto& Xf = X.ViewComponent("face", true);
  const auto& Xc = X.ViewComponent("cell");

  {
    auto Yf = Y.ViewComponent("face", true);
    auto Yc = Y.ViewComponent("cell");

    const auto& map = Y.getMap();

    for (int c = 0; c != ncells_owned; ++c) {
      Kokkos::View<AmanziMesh::Entity_ID*> faces;
      mesh_->cell_get_faces(c, faces);
      int nfaces = faces.size();

      int npoints(0);
      // TODO correct for map (blockmap)
      for (int n = 0; n != nfaces; ++n) npoints += 0;//map.ElementSize(faces[n]);

      int m(0);
      WhetStone::DenseVector v(npoints + 1), av(npoints + 1);

      for (int n = 0; n != nfaces; ++n) {
        int f = faces[n];
        // TODO not working for tpetra map (blockmap yes)
        int first = 0;//map.FirstPointInElement(f);
        for (int k = 0; k < 0; ++k){//map.ElementSize(f); ++k) {
          v(m++) = Xf(0,first + k);
        }
      }
      v(npoints) = Xc(0,c);

      // TODO build a matrix based on view<double**>
      //const WhetStone::DenseMatrix& Acell = op.matrices[c];
      //Acell.elementWiseMultiply(v, av, false);

      m = 0;
      for (int n = 0; n != nfaces; ++n) {
        int f = faces[n];
        // TODO not working for tpetra map (blockmap yes)
        int first = 0;//map.FirstPointInElement(f);
        for (int k = 0; k < 0; ++k){//map.ElementSize(f); ++k) {
          Yf(0,first + k) += av(m++);
        }
      }
      Yc(0,c) += av(npoints);
    }
  }

  //Y.GatherGhostedToMaster(Add);
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
  AMANZI_ASSERT(op.data.extent(0) == ncells_owned);

  Y.putScalarGhosted(0.);
  X.ScatterMasterToGhosted();
  const auto& Xf = X.ViewComponent("face", true);

  {
    auto Yf = Y.ViewComponent("face", true);

    for (int c = 0; c != ncells_owned; ++c) {
      Kokkos::View<AmanziMesh::Entity_ID*> faces;
      mesh_->cell_get_faces(c, faces);
      int nfaces = faces.size();

      WhetStone::DenseVector v(nfaces), av(nfaces);
      for (int n = 0; n != nfaces; ++n) { v(n) = Xf(0,faces[n]); }

      // TODO  Build matrix based on view<double**>
      //const WhetStone::DenseMatrix& Acell = op.matrices[c];
      //Acell.elementWiseMultiply(v, av, false);

      for (int n = 0; n != nfaces; ++n) { Yf(0,faces[n]) += av(n); }
    }
  }
  // TODO 
  //Y.GatherGhostedToMaster(Add);
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
  int nsurf_cells = op.surf_mesh->num_entities(
    AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  // \TODO extract diag from op 
  //AMANZI_ASSERT(op.diag->getLocalLength() == nsurf_cells);

  const auto& Xf = X.ViewComponent("face", false);
  auto Yf = Y.ViewComponent("face", false);
  for (int sc = 0; sc != nsurf_cells; ++sc) {
    int f = op.surf_mesh->entity_get_parent(AmanziMesh::CELL, sc);
    // TODO extract diag from op
    //Yf(0,f) += (*op.diag)(0,sc) * Xf(0,f);
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
  int nsurf_faces = op.surf_mesh->num_entities(
    AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  AMANZI_ASSERT(op.data.extent(0) == nsurf_faces);

  X.ScatterMasterToGhosted();
  const auto& Xf = X.ViewComponent("face", true);

  Y.putScalarGhosted(0.);
  auto Yf = Y.ViewComponent("face", true);

  for (int sf = 0; sf != nsurf_faces; ++sf) {
    Kokkos::View<AmanziMesh::Entity_ID*> cells;
    op.surf_mesh->face_get_cells(sf, AmanziMesh::Parallel_type::ALL, cells);
    int ncells = cells.extent(0);

    WhetStone::DenseVector v(ncells), av(ncells);
    for (int n = 0; n != ncells; ++n) {
      v(n) = Xf(0,op.surf_mesh->entity_get_parent(AmanziMesh::CELL, cells[n]));
    }

    // TODO extract matrix from view<double**>
    //const WhetStone::DenseMatrix& Aface = op.matrices[sf];
    //Aface.elementWiseMultiply(v, av, false);

    for (int n = 0; n != ncells; ++n) {
      Yf(0,op.surf_mesh->entity_get_parent(AmanziMesh::CELL, cells[n])) +=
        av(n);
    }
  }
  Y.GatherGhostedToMaster("face");
  return 0;
}


/* ******************************************************************
 * Visit methods for symbolic assemble: FaceCell
 ****************************************************************** */
void
Operator_FaceCell::SymbolicAssembleMatrixOp(const Op_Cell_FaceCell& op,
                                            const SuperMap& map, GraphFE& graph,
                                            int my_block_row,
                                            int my_block_col) const
{
  std::vector<int> lid_r(2 * cell_max_faces + 1);
  std::vector<int> lid_c(2 * cell_max_faces + 1);

  // ELEMENT: cell, DOFS: cell and face
  const auto face_row_inds =
    map.GhostIndices(my_block_row, "face", 0);
  const auto face_col_inds =
    map.GhostIndices(my_block_col, "face", 0);
  const auto cell_row_inds =
    map.GhostIndices(my_block_row, "cell", 0);
  const auto cell_col_inds =
    map.GhostIndices(my_block_col, "cell", 0);

  auto face_gh_map =
    map.ComponentGhostedMap(my_block_row, "face");
  auto cell_map =
    map.ComponentMap(my_block_row, "cell");

  int ierr(0);
  for (int c = 0; c != ncells_owned; ++c) {
    Kokkos::View<AmanziMesh::Entity_ID*> faces;
    mesh_->cell_get_faces(c, faces);
    int nfaces = faces.size();

    int k = 0;
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      // TODO Not working for tpetra map (blockmap yes)
      int face_dof_size = 0;//face_gh_map->ElementSize(f);
      int first = 0;//face_gh_map->FirstPointInElement(f);

      for (int m = 0; m != face_dof_size; ++m) {
        lid_r[k] = face_row_inds[first + m];
        lid_c[k] = face_col_inds[first + m];
        k++;
      }
    }

    // TODO Not working for tpetra map (blockmap yes)
    int cell_dof_size = 0;//cell_map->ElementSize(c);
    int first = 0;//cell_map->FirstPointInElement(c);
    for (int m = 0; m != cell_dof_size; ++m) {
      lid_r[k] = cell_row_inds[first + m];
      lid_c[k] = cell_col_inds[first + m];
      k++;
    }
    // TODO 
    //ierr |= graph.InsertMyIndices(k, lid_r.data(), k, lid_c.data());
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
 * Visit methods for symbolic assemble: Face
 ****************************************************************** */
void
Operator_FaceCell::SymbolicAssembleMatrixOp(const Op_Cell_Face& op,
                                            const SuperMap& map, GraphFE& graph,
                                            int my_block_row,
                                            int my_block_col) const
{
  std::vector<int> lid_r(2 * cell_max_faces);
  std::vector<int> lid_c(2 * cell_max_faces);

  // ELEMENT: cell, DOFS: face
  const auto face_row_inds =
    map.GhostIndices(my_block_row, "face", 0);
  const auto face_col_inds =
    map.GhostIndices(my_block_col, "face", 0);

  auto face_gh_map =
    map.ComponentGhostedMap(my_block_row, "face");

  int ierr(0);
  for (int c = 0; c != ncells_owned; ++c) {
    Kokkos::View<AmanziMesh::Entity_ID*> faces;
    mesh_->cell_get_faces(c, faces);
    int nfaces = faces.size();

    int k = 0;
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      // TODO Do not work with Tpetra map (blockmap yes)
      int face_dof_size = 0;//face_gh_map->ElementSize(f);
      int first = 0;//face_gh_map->FirstPointInElement(f);

      for (int m = 0; m != face_dof_size; ++m) {
        lid_r[k] = face_row_inds[first + m];
        lid_c[k] = face_col_inds[first + m];
        k++;
      }
    }
    // TODO 
    //ierr |= graph.InsertMyIndices(k, lid_r.data(), k, lid_c.data());
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
 * Visit methods for symbolic assemble: Surface
 ****************************************************************** */
void
Operator_FaceCell::SymbolicAssembleMatrixOp(
  const Op_SurfaceCell_SurfaceCell& op, const SuperMap& map, GraphFE& graph,
  int my_block_row, int my_block_col) const
{
  int nsurf_cells = op.surf_mesh->num_entities(
    AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  // ELEMENT: cell, DOFS: cell and face
  const auto face_row_inds =
    map.GhostIndices(my_block_row, "face", 0);
  const auto face_col_inds =
    map.GhostIndices(my_block_col, "face", 0);

  int ierr = 0;
  for (int sc = 0; sc != nsurf_cells; ++sc) {
    int lid_r =
      face_row_inds[op.surf_mesh->entity_get_parent(AmanziMesh::CELL, sc)];
    int lid_c =
      face_col_inds[op.surf_mesh->entity_get_parent(AmanziMesh::CELL, sc)];
    // TODO 
    //ierr |= graph.InsertMyIndices(lid_r, 1, &lid_c);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
 * Visit methods for symbolic assemble: Surface
 ****************************************************************** */
void
Operator_FaceCell::SymbolicAssembleMatrixOp(
  const Op_SurfaceFace_SurfaceCell& op, const SuperMap& map, GraphFE& graph,
  int my_block_row, int my_block_col) const
{
  int nsurf_faces = op.surf_mesh->num_entities(
    AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int lid_r[2];
  int lid_c[2];

  // ELEMENT: cell, DOFS: cell and face
  const auto face_row_inds =
    map.GhostIndices(my_block_row, "face", 0);
  const auto face_col_inds =
    map.GhostIndices(my_block_col, "face", 0);

  int ierr = 0;
  for (int sf = 0; sf != nsurf_faces; ++sf) 
  {
    Kokkos::View<AmanziMesh::Entity_ID*> cells;
    op.surf_mesh->face_get_cells(sf, AmanziMesh::Parallel_type::ALL, cells);
    int ncells = cells.extent(0);
    for (int n = 0; n != ncells; ++n) {
      int f = op.surf_mesh->entity_get_parent(AmanziMesh::CELL, cells[n]);
      lid_r[n] = face_row_inds[f];
      lid_c[n] = face_col_inds[f];
    }
    // TODO 
    //ierr |= graph.InsertMyIndices(ncells, lid_r, ncells, lid_c);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
 * Visit methods for symbolic assemble: Coupling
 ****************************************************************** */
void
Operator_FaceCell::SymbolicAssembleMatrixOp(const Op_Diagonal& op,
                                            const SuperMap& map, GraphFE& graph,
                                            int my_block_row,
                                            int my_block_col) const
{
  std::string row_name = op.row_compname();
  std::string col_name = op.col_compname();

  const auto row_gids =
    map.GhostIndices(my_block_row, op.row_compname(), 0);
  const auto col_gids =
    map.GhostIndices(my_block_col, op.col_compname(), 0);

  const auto& col_lids = op.col_inds();
  const auto& row_lids = op.row_inds();

  std::vector<int> lid_r, lid_c;

  int ierr(0);
  for (int n = 0; n != col_lids.size(); ++n) {
    int ndofs = col_lids[n].size();

    lid_r.clear();
    lid_c.clear();

    for (int i = 0; i != ndofs; ++i) {
      lid_r.push_back(row_gids[row_lids[n][i]]);
      lid_c.push_back(col_gids[col_lids[n][i]]);
    }
    // TODO
    //ierr |= graph.InsertMyIndices(ndofs, lid_r.data(), ndofs, lid_c.data());
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
 * Visit methods for assemble: FaceCell
 ****************************************************************** */
void
Operator_FaceCell::AssembleMatrixOp(const Op_Cell_FaceCell& op,
                                    const SuperMap& map, MatrixFE& mat,
                                    int my_block_row, int my_block_col) const
{
  AMANZI_ASSERT(op.data.extent(0) == ncells_owned);

  std::vector<int> lid_r(2 * cell_max_faces + 1);
  std::vector<int> lid_c(2 * cell_max_faces + 1);

  // ELEMENT: cell, DOFS: face and cell
  const auto face_row_inds =
    map.GhostIndices(my_block_row, "face", 0);
  const auto face_col_inds =
    map.GhostIndices(my_block_col, "face", 0);
  const auto cell_row_inds =
    map.GhostIndices(my_block_row, "cell", 0);
  const auto cell_col_inds =
    map.GhostIndices(my_block_col, "cell", 0);

  auto face_gh_map =
    map.ComponentGhostedMap(my_block_row, "face");
  auto cell_map =
    map.ComponentMap(my_block_row, "cell");

  int ierr(0);
  for (int c = 0; c != ncells_owned; ++c) {
    Kokkos::View<AmanziMesh::Entity_ID*> faces;
    mesh_->cell_get_faces(c, faces);

    int nfaces = faces.size();
    int k = 0;
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      // TODO not working for tpetra map (blockmap yes)
      int face_dof_size = 0;//face_gh_map->ElementSize(f);
      int first = 0;//face_gh_map->FirstPointInElement(f);

      for (int m = 0; m != face_dof_size; ++m) {
        lid_r[k] = face_row_inds[first + m];
        lid_c[k] = face_col_inds[first + m];
        k++;
      }
    }

    // TODO not working for tpetra map (blockmap yes)
    int cell_dof_size = 0;//cell_map->ElementSize(c);
    int first = 0;//cell_map->FirstPointInElement(c);
    for (int m = 0; m != cell_dof_size; ++m) {
      lid_r[k] = cell_row_inds[first + m];
      lid_c[k] = cell_col_inds[first + m];
      k++;
    }
    // TODO 
    //ierr |= mat.SumIntoMyValues(lid_r.data(), lid_c.data(), op.matrices[c]);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
 * Visit methods for assemble: Face
 ****************************************************************** */
void
Operator_FaceCell::AssembleMatrixOp(const Op_Cell_Face& op, const SuperMap& map,
                                    MatrixFE& mat, int my_block_row,
                                    int my_block_col) const
{
  AMANZI_ASSERT(op.data.extent(0) == ncells_owned);

  std::vector<int> lid_r(cell_max_faces + 1);
  std::vector<int> lid_c(cell_max_faces + 1);

  // ELEMENT: cell, DOFS: face and cell
  const auto face_row_inds =
    map.GhostIndices(my_block_row, "face", 0);
  const auto face_col_inds =
    map.GhostIndices(my_block_col, "face", 0);

  auto face_gh_map =
    map.ComponentGhostedMap(my_block_row, "face");

  int ierr(0);
  for (int c = 0; c != ncells_owned; ++c) {
    Kokkos::View<AmanziMesh::Entity_ID*> faces;
    mesh_->cell_get_faces(c, faces);
    int nfaces = faces.extent(0);

    int k = 0;
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      // TODO not available for map (blockmap yes)
      int face_dof_size = 0;// face_gh_map->ElementSize(f);
      int first = 0;//face_gh_map->FirstPointInElement(f);

      for (int m = 0; m != face_dof_size; ++m) {
        lid_r[k] = face_row_inds[first + m];
        lid_c[k] = face_col_inds[first + m];
        k++;
      }
    }
    // TODO 
    //ierr |= mat.SumIntoMyValues(lid_r.data(), lid_c.data(), op.matrices[c]);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
 * Visit methods for assemble: Surface
 ****************************************************************** */
void
Operator_FaceCell::AssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op,
                                    const SuperMap& map, MatrixFE& mat,
                                    int my_block_row, int my_block_col) const
{
  int nsurf_cells = op.surf_mesh->num_entities(
    AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  // ELEMENT: cell, DOFS: cell and face
  const auto face_row_inds =
    map.GhostIndices(my_block_row, "face", 0);
  const auto face_col_inds =
    map.GhostIndices(my_block_col, "face", 0);

  int ierr = 0;
  for (int sc = 0; sc != nsurf_cells; ++sc) {
    int lid_r =
      face_row_inds[op.surf_mesh->entity_get_parent(AmanziMesh::CELL, sc)];
    int lid_c =
      face_col_inds[op.surf_mesh->entity_get_parent(AmanziMesh::CELL, sc)];
    // TODO 
    //ierr |= mat.SumIntoMyValues(lid_r, 1, &(*op.diag)[0][sc], &lid_c);
  }
  AMANZI_ASSERT(!ierr);
}


void
Operator_FaceCell::AssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op,
                                    const SuperMap& map, MatrixFE& mat,
                                    int my_block_row, int my_block_col) const
{
  int nsurf_faces = op.surf_mesh->num_entities(
    AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int lid_r[2];
  int lid_c[2];

  // ELEMENT: cell, DOFS: cell and face
  const auto face_row_inds =
    map.GhostIndices(my_block_row, "face", 0);
  const auto face_col_inds =
    map.GhostIndices(my_block_col, "face", 0);

  int ierr = 0;
  for (int sf = 0; sf != nsurf_faces; ++sf) {
    Kokkos::View<AmanziMesh::Entity_ID*> cells;
    op.surf_mesh->face_get_cells(sf, AmanziMesh::Parallel_type::ALL, cells);
    int ncells = cells.extent(0);
    for (int n = 0; n != ncells; ++n) {
      int f = op.surf_mesh->entity_get_parent(AmanziMesh::CELL, cells[n]);
      lid_r[n] = face_row_inds[f];
      lid_c[n] = face_col_inds[f];
    }
    // TODO 
    //ierr |= mat.SumIntoMyValues(lid_r, lid_c, op.matrices[sf]);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
 * Visit methods for assemble: Coupling
 ****************************************************************** */
void
Operator_FaceCell::AssembleMatrixOp(const Op_Diagonal& op, const SuperMap& map,
                                    MatrixFE& mat, int my_block_row,
                                    int my_block_col) const
{
  const auto row_gids = 
    map.GhostIndices(my_block_row, op.row_compname(), 0);
  const auto col_gids = 
    map.GhostIndices(my_block_col, op.col_compname(), 0);

  const auto& col_lids = op.col_inds();
  const auto& row_lids = op.row_inds();

  std::vector<int> lid_r, lid_c;

  int ierr(0);
  for (int n = 0; n != col_lids.size(); ++n) {
    int ndofs = col_lids[n].size();

    lid_r.clear();
    lid_c.clear();

    for (int i = 0; i != ndofs; ++i) {
      lid_r.push_back(row_gids[row_lids[n][i]]);
      lid_c.push_back(col_gids[col_lids[n][i]]);
    }
    // TODO 
    //ierr |= mat.SumIntoMyValues(lid_r.data(), lid_c.data(), op.matrices[n]);
  }
  AMANZI_ASSERT(!ierr);
}

} // namespace Operators
} // namespace Amanzi
