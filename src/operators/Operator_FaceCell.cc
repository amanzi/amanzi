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
int Operator_FaceCell::ApplyMatrixFreeOp(const Op_Cell_FaceCell& op,
                                         const CompositeVector& X, CompositeVector& Y) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

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
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

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
  int nsurf_cells = op.surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  AMANZI_ASSERT(op.diag->MyLength() == nsurf_cells);

  const Epetra_MultiVector& Xf = *X.ViewComponent("face", false);
  Epetra_MultiVector& Yf = *Y.ViewComponent("face", false);
  for (int sc = 0; sc != nsurf_cells; ++sc) {
    int f = op.surf_mesh->entity_get_parent(AmanziMesh::CELL, sc);
    Yf[0][f] += (*op.diag)[0][sc] * Xf[0][f];
  } 
  return 0;
}


/* ******************************************************************
* visit method for apply surface cells into subsurface faces
****************************************************************** */
int Operator_FaceCell::ApplyMatrixFreeOp(const Op_SurfaceFace_SurfaceCell& op,
                                         const CompositeVector& X, CompositeVector& Y) const
{
  int nsurf_faces = op.surf_mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  AMANZI_ASSERT(op.matrices.size() == nsurf_faces);

  X.ScatterMasterToGhosted();
  const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

  Y.PutScalarGhosted(0.);
  Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);
  
  AmanziMesh::Entity_ID_List cells;
  for (int sf = 0; sf != nsurf_faces; ++sf) {
    op.surf_mesh->face_get_cells(sf, AmanziMesh::Parallel_type::ALL, &cells);
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
* Visit methods for Apply with variable number of DOFs (aka points):
* apply the local matrices directly as schemas match.
****************************************************************** */
int Operator_FaceCell::ApplyMatrixFreeOpVariableDOFs(
    const Op_Cell_FaceCell& op,
    const CompositeVector& X, CompositeVector& Y) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

  Y.PutScalarGhosted(0.0);
  X.ScatterMasterToGhosted();
  const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);
  const Epetra_MultiVector& Xc = *X.ViewComponent("cell");

  {
    Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);
    Epetra_MultiVector& Yc = *Y.ViewComponent("cell");

    const auto& map = Yf.Map();

    AmanziMesh::Entity_ID_List faces;
    for (int c = 0; c != ncells_owned; ++c) {
      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      int npoints(0);
      for (int n = 0; n != nfaces; ++n)
        npoints += map.ElementSize(faces[n]);

      int m(0);
      WhetStone::DenseVector v(npoints + 1), av(npoints + 1);

      for (int n = 0; n != nfaces; ++n) {
        int f = faces[n];
        int first = map.FirstPointInElement(f);
        for (int k = 0; k < map.ElementSize(f); ++k) {
          v(m++) = Xf[0][first + k];
        }
      }
      v(npoints) = Xc[0][c];

      const WhetStone::DenseMatrix& Acell = op.matrices[c];
      Acell.Multiply(v, av, false);

      m = 0;
      for (int n = 0; n != nfaces; ++n) {
        int f = faces[n];
        int first = map.FirstPointInElement(f);
        for (int k = 0; k < map.ElementSize(f); ++k) {
          Yf[0][first + k] += av(m++);
        }
      }
      Yc[0][c] += av(npoints);
    } 
  }

  Y.GatherGhostedToMaster(Add);
  return 0;
}


/* ******************************************************************
* Visit methods for symbolic assemble: FaceCell
****************************************************************** */
void Operator_FaceCell::SymbolicAssembleMatrixOp(const Op_Cell_FaceCell& op,
                                                 const SuperMap& map, GraphFE& graph,
                                                 int my_block_row, int my_block_col, bool multi_domain) const
{
  std::vector<int> lid_r(2*cell_max_faces + 1);
  std::vector<int> lid_c(2*cell_max_faces + 1);

  std::string face_name = "face";
  std::string cell_name = "cell";
  int row_pos = my_block_row;
  int col_pos = my_block_col;

  if (multi_domain) {
    face_name = face_name + "-" + std::to_string(my_block_row);
    cell_name = cell_name + "-" + std::to_string(my_block_col);
    row_pos = 0;
    col_pos = 0;
  }
    
  // ELEMENT: cell, DOFS: cell and face
  Teuchos::RCP<const Epetra_BlockMap> face_gh_map = map.BaseGhostedMap(face_name);
  Teuchos::RCP<const Epetra_BlockMap> cell_gh_map = map.BaseGhostedMap(cell_name);
  Teuchos::RCP<const Epetra_BlockMap> face_map = map.BaseMap(face_name);
  int nface_points_owned = face_map->NumMyPoints();
  int num_dof_faces(0), num_dof_cells(0);
  num_dof_faces = map.NumDofs(face_name);
  num_dof_cells = map.NumDofs(cell_name);

  int face_row_offset = map.Offset(face_name);
  int face_gh_offset = map.GhostedOffset(face_name);
  int cell_row_offset = map.Offset(cell_name);
  
  int ierr(0);
  AmanziMesh::Entity_ID_List faces;
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    int k = 0;
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      int face_dof_size = face_gh_map->ElementSize(f);
      int first = face_gh_map->FirstPointInElement(f);

      first = (faces[n] < nfaces_owned) ? first : (first - nface_points_owned);
      int offset = (faces[n] < nfaces_owned) ? face_row_offset : face_gh_offset;
      
      for (int m = 0; m != face_dof_size; ++m) {
        lid_r[k] = offset + (first + m) * num_dof_faces + row_pos;
        lid_c[k] = offset + (first + m) * num_dof_faces + col_pos;
        k++;
      }
    }

    int cell_dof_size = cell_gh_map->ElementSize(c);
    int first = cell_gh_map->FirstPointInElement(c);
    for (int m = 0; m != cell_dof_size; ++m) {
      lid_r[k] = cell_row_offset + (first + m) * num_dof_cells + row_pos;
      lid_c[k] = cell_row_offset + (first + m) * num_dof_cells + col_pos;
      k++;
    }
    
    ierr |= graph.InsertMyIndices(k, lid_r.data(), k, lid_c.data());
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for symbolic assemble: Face
****************************************************************** */
void Operator_FaceCell::SymbolicAssembleMatrixOp(const Op_Cell_Face& op,
                                                 const SuperMap& map, GraphFE& graph,
                                                 int my_block_row, int my_block_col, bool multi_domain) const
{
  std::vector<int> lid_r(2*cell_max_faces);
  std::vector<int> lid_c(2*cell_max_faces);

  // ELEMENT: cell, DOFS: face
  const std::vector<int>& face_row_inds = map.GhostIndices("face", my_block_row);
  const std::vector<int>& face_col_inds = map.GhostIndices("face", my_block_col);

  Teuchos::RCP<const Epetra_BlockMap> face_gh_map = map.BaseGhostedMap("face");
  Teuchos::RCP<const Epetra_BlockMap> face_map = map.BaseMap("face");
  int nface_points_owned = face_map->NumMyPoints();
  int num_dof_faces = map.NumDofs("face");
  
  int face_row_offset = map.Offset("face");
  int face_gh_offset = map.GhostedOffset("face");

  int ierr(0);
  AmanziMesh::Entity_ID_List faces;
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    int k = 0;
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      int face_dof_size = face_map->ElementSize(f);
      int first = face_gh_map->FirstPointInElement(f);

      first = (faces[n] < nfaces_owned) ? first : (first - nface_points_owned);
      int offset = (faces[n] < nfaces_owned) ? face_row_offset : face_gh_offset;
      
      for (int m = 0; m != face_dof_size; ++m) {
        lid_r[k] = offset + (first + m) * num_dof_faces + my_block_row;
        lid_c[k] = offset + (first + m) * num_dof_faces + my_block_col;
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
void Operator_FaceCell::SymbolicAssembleMatrixOp(
    const Op_SurfaceCell_SurfaceCell& op,
    const SuperMap& map, GraphFE& graph,
    int my_block_row, int my_block_col, bool multi_domain) const
{
  int nsurf_cells = op.surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices("face", my_block_row);
  const std::vector<int>& face_col_inds = map.GhostIndices("face", my_block_col);

  int ierr = 0;
  for (int sc = 0; sc != nsurf_cells; ++sc) {
    int lid_r = face_row_inds[op.surf_mesh->entity_get_parent(AmanziMesh::CELL, sc)];
    int lid_c = face_col_inds[op.surf_mesh->entity_get_parent(AmanziMesh::CELL, sc)];
    ierr |= graph.InsertMyIndices(lid_r, 1, &lid_c);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for symbolic assemble: Surface
****************************************************************** */
void Operator_FaceCell::SymbolicAssembleMatrixOp(
    const Op_SurfaceFace_SurfaceCell& op,
    const SuperMap& map, GraphFE& graph,
    int my_block_row, int my_block_col, bool multi_domain) const
{
  int nsurf_faces = op.surf_mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int lid_r[2];
  int lid_c[2];

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices("face", my_block_row);
  const std::vector<int>& face_col_inds = map.GhostIndices("face", my_block_col);

  int ierr = 0;
  AmanziMesh::Entity_ID_List cells;
  for (int sf = 0; sf != nsurf_faces; ++sf) {
    op.surf_mesh->face_get_cells(sf, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();
    for (int n = 0; n != ncells; ++n) {
      int f = op.surf_mesh->entity_get_parent(AmanziMesh::CELL, cells[n]);
      lid_r[n] = face_row_inds[f];
      lid_c[n] = face_col_inds[f];
    }

    ierr |= graph.InsertMyIndices(ncells, lid_r, ncells, lid_c);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for symbolic assemble: Coupling
****************************************************************** */
void Operator_FaceCell::SymbolicAssembleMatrixOp(const Op_Diagonal& op,
                                                 const SuperMap& map, GraphFE& graph,
                                                 int my_block_row, int my_block_col, bool multi_domain) const
{
  std::string row_name = op.row_compname();
  std::string col_name = op.col_compname();
  int row_pos = my_block_row;
  int col_pos = my_block_col;

  if (multi_domain) {
    row_name = row_name + "-" + std::to_string(my_block_row);
    col_name = col_name + "-" + std::to_string(my_block_col);
    row_pos = 0;
    col_pos = 0;
  }
        
  const std::vector<int>& row_gids = map.GhostIndices(row_name, row_pos);
  const std::vector<int>& col_gids = map.GhostIndices(col_name, col_pos);

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
    ierr |= graph.InsertMyIndices(ndofs, lid_r.data(), ndofs, lid_c.data());
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble: FaceCell
****************************************************************** */
void Operator_FaceCell::AssembleMatrixOp(const Op_Cell_FaceCell& op,
                                         const SuperMap& map, MatrixFE& mat,
                                         int my_block_row, int my_block_col, bool multi_domain) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

  std::vector<int> lid_r(2*cell_max_faces + 1);
  std::vector<int> lid_c(2*cell_max_faces + 1);


  std::string face_name = "face";
  std::string cell_name = "cell";
  int row_pos = my_block_row;
  int col_pos = my_block_col;

  if (multi_domain) {
    face_name = face_name + "-" + std::to_string(my_block_row);
    cell_name = cell_name + "-" + std::to_string(my_block_col);
    row_pos = 0;
    col_pos = 0;
  }

  // ELEMENT: cell, DOFS: face and cell
  const std::vector<int>& face_row_inds = map.GhostIndices(face_name, row_pos);
  const std::vector<int>& face_col_inds = map.GhostIndices(face_name, col_pos);
  const std::vector<int>& cell_row_inds = map.GhostIndices(cell_name, row_pos);
  const std::vector<int>& cell_col_inds = map.GhostIndices(cell_name, col_pos);        

  Teuchos::RCP<const Epetra_BlockMap> face_gh_map = map.BaseGhostedMap(face_name);
  Teuchos::RCP<const Epetra_BlockMap> cell_gh_map = map.BaseGhostedMap(cell_name);
  Teuchos::RCP<const Epetra_BlockMap> face_map = map.BaseMap(face_name);

  int nface_points_owned = face_map->NumMyPoints();
  int num_dof_faces = map.NumDofs(face_name);
  int num_dof_cells = map.NumDofs(cell_name);
  int face_row_offset = map.Offset(face_name);
  int face_gh_offset = map.GhostedOffset(face_name);
  int cell_row_offset = map.Offset(cell_name);

  int ierr(0);
  AmanziMesh::Entity_ID_List faces;
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    
    int nfaces = faces.size();
    int k = 0;
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      int face_dof_size = face_gh_map->ElementSize(f);
      int first = face_gh_map->FirstPointInElement(f);

      first = (faces[n] < nfaces_owned) ? first : (first - nface_points_owned);
      int offset = (faces[n] < nfaces_owned) ? face_row_offset : face_gh_offset;
      
      for (int m = 0; m != face_dof_size; ++m) {
        lid_r[k] = offset + (first + m) * num_dof_faces + row_pos;
        lid_c[k] = offset + (first + m) * num_dof_faces + col_pos;
        k++;
      }      
    }

    int cell_dof_size = cell_gh_map->ElementSize(c);
    int first = cell_gh_map->FirstPointInElement(c);
    for (int m = 0; m != cell_dof_size; ++m) {
      lid_r[k] = cell_row_offset + (first + m) * num_dof_cells + row_pos;
      lid_c[k] = cell_row_offset + (first + m) * num_dof_cells + col_pos;
      k++;
    }
    
    ierr |= mat.SumIntoMyValues(lid_r.data(), lid_c.data(), op.matrices[c]);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble: Face
****************************************************************** */
void Operator_FaceCell::AssembleMatrixOp(const Op_Cell_Face& op,
                                         const SuperMap& map, MatrixFE& mat,
                                         int my_block_row, int my_block_col, bool multi_domain) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

  std::vector<int> lid_r(cell_max_faces + 1);
  std::vector<int> lid_c(cell_max_faces + 1);

  // ELEMENT: cell, DOFS: face and cell
  Teuchos::RCP<const Epetra_BlockMap> face_gh_map = map.BaseGhostedMap("face");
  Teuchos::RCP<const Epetra_BlockMap> face_map = map.BaseMap("face");
  int nface_points_owned = face_map->NumMyPoints();
  int num_dof_faces = map.NumDofs("face");
  
  int face_row_offset = map.Offset("face");
  int face_gh_offset = map.GhostedOffset("face");
  
  int ierr(0);
  AmanziMesh::Entity_ID_List faces;
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_faces(c, &faces);
    
    int nfaces = faces.size();

    int k = 0;
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      int face_dof_size = face_map->ElementSize(f);
      int first = face_gh_map->FirstPointInElement(f);

      first = (faces[n] < nfaces_owned) ? first : (first - nface_points_owned);
      int offset = (faces[n] < nfaces_owned) ? face_row_offset : face_gh_offset;
      
      for (int m = 0; m != face_dof_size; ++m) {
        lid_r[k] = offset + (first + m) * num_dof_faces + my_block_row;
        lid_c[k] = offset + (first + m) * num_dof_faces + my_block_col;
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
void Operator_FaceCell::AssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op,
                                         const SuperMap& map, MatrixFE& mat,
                                         int my_block_row, int my_block_col, bool multi_domain) const
{
  int nsurf_cells = op.surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices("face", my_block_row);
  const std::vector<int>& face_col_inds = map.GhostIndices("face", my_block_col);

  int ierr = 0;
  for (int sc = 0; sc != nsurf_cells; ++sc) {
    int lid_r = face_row_inds[op.surf_mesh->entity_get_parent(AmanziMesh::CELL, sc)];
    int lid_c = face_col_inds[op.surf_mesh->entity_get_parent(AmanziMesh::CELL, sc)];
    ierr |= mat.SumIntoMyValues(lid_r, 1, &(*op.diag)[0][sc], &lid_c);
  }
  AMANZI_ASSERT(!ierr);
}


void Operator_FaceCell::AssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op,
                                         const SuperMap& map, MatrixFE& mat,
                                         int my_block_row, int my_block_col, bool multi_domain) const
{
  int nsurf_faces = op.surf_mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int lid_r[2];
  int lid_c[2];

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices("face", my_block_row);
  const std::vector<int>& face_col_inds = map.GhostIndices("face", my_block_col);

  int ierr = 0;
  AmanziMesh::Entity_ID_List cells;
  for (int sf = 0; sf != nsurf_faces; ++sf) {
    op.surf_mesh->face_get_cells(sf, AmanziMesh::Parallel_type::ALL, &cells);
    int ncells = cells.size();
    for (int n = 0; n != ncells; ++n) {
      int f = op.surf_mesh->entity_get_parent(AmanziMesh::CELL, cells[n]);
      lid_r[n] = face_row_inds[f];
      lid_c[n] = face_col_inds[f];
    }

    ierr |= mat.SumIntoMyValues(lid_r, lid_c, op.matrices[sf]);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble: Coupling
****************************************************************** */
void Operator_FaceCell::AssembleMatrixOp(const Op_Diagonal& op,
                                         const SuperMap& map, MatrixFE& mat,
                                         int my_block_row, int my_block_col, bool multi_domain) const
{
  std::string row_name = op.row_compname();
  std::string col_name = op.col_compname();
  int row_pos = my_block_row;
  int col_pos = my_block_col;

  if (multi_domain) {
    row_name = row_name + "-" + std::to_string(my_block_row);
    col_name = col_name + "-" + std::to_string(my_block_col);
    row_pos = 0;
    col_pos = 0;    
  }
  
  const std::vector<int>& row_gids = map.GhostIndices(row_name, row_pos);
  const std::vector<int>& col_gids = map.GhostIndices(col_name, col_pos);

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
      // std::cout << row_lids[n][i]<<" "<<col_lids[n][i]<<" : "<<
      //   row_gids[row_lids[n][i]] << " "<<col_gids[col_lids[n][i]]<<"\n";      
    }

    ierr |= mat.SumIntoMyValues(lid_r.data(), lid_c.data(), op.matrices[n]);
  }
  AMANZI_ASSERT(!ierr);
}

}  // namespace Operators
}  // namespace Amanzi



