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
//#include "Op_Diagonal.hh"
//#include "Op_SurfaceCell_SurfaceCell.hh"
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
  std::cout<<"Operator_FaceCell::ApplyMatrixFreeOp"<<std::endl;
  std::cout<<"nfaces_owned: "<<nfaces_owned<<" op.A.size(): "<<op.A.size()<<std::endl;
  AMANZI_ASSERT(op.A.size() == ncells_owned);
  auto Xf = X.ViewComponent("face", true);
  auto Xc = X.ViewComponent("cell", true);

  auto Yf = Y.ViewComponent("face", true);
  auto Yc = Y.ViewComponent("cell", true);
  
  const auto& map = Y.getMap();
  const AmanziMesh::Mesh* mesh = mesh_.get();

  auto local_A = op.A; 
  auto local_Av = op.Av; 
  auto local_v = op.v; 
  // Allocate the first time 
  if (local_v.size() != local_A.size()) {
    op.PreallocateWorkVectors();
    local_v = op.v;
    local_Av = op.Av;
  }

  Kokkos::parallel_for(
    "Operator_FaceCell::ApplyMatrixFreeOp",
    ncells_owned, 
    KOKKOS_LAMBDA(const int& c){
      AmanziMesh::Entity_ID_View faces;
      mesh->cell_get_faces(c, faces);
      int nfaces = faces.size();

      //int npoints(0);
      //for (int n = 0; n != nfaces; ++n)
      //  npoints += map.ElementSize(faces[n]);

      WhetStone::DenseVector<Amanzi::DeviceOnlyMemorySpace> lv = getFromCSR<WhetStone::DenseVector>(local_v,c);
      WhetStone::DenseVector<Amanzi::DeviceOnlyMemorySpace> lAv = getFromCSR<WhetStone::DenseVector>(local_Av,c);
      WhetStone::DenseMatrix<Amanzi::DeviceOnlyMemorySpace> lA = getFromCSR<WhetStone::DenseMatrix>(local_A,c);
      //WhetStone::DenseVector<DeviceOnlyMemorySpace> v(npoints + 1), av(npoints + 1);

      int m(0);
      for (int n = 0; n != nfaces; ++n) {
      //  int f = faces[n];
      //  assert(false); 
        //int first = map.FirstPointInElement(f);
        //for (int k = 0; k < map.ElementSize(f); ++k) {
        lv(m++) = Xf(0,n);
        //}
      }
      lv(nfaces) = Xc(0,c); 
      //lv(npoints) = Xc(0,c);

      lA.Multiply(lv, lAv, false);

      m = 0;
      for (int n = 0; n != nfaces; ++n) {
        int f = faces[n];
        Kokkos::atomic_add(&Yf(0,n), lAv(m++));
      }
      Yc(0,c) += lAv(nfaces);
    });  
  return 0;
}

#if 0 
/* ******************************************************************
* Apply the local matrices directly as schemas match.
****************************************************************** */
int Operator_FaceCell::ApplyMatrixFreeOp(const Op_Cell_Face& op,
                                         const CompositeVector& X, CompositeVector& Y) const
{

  auto Xf = X.ViewComponent("face", true);
  auto Yf = Y.ViewComponent("face", true);

  const AmanziMesh::Mesh* mesh = mesh_.get();

  auto local_A = op.A; 
  auto local_Av = op.Av; 
  auto local_v = op.v; 
  // Allocate the first time 
  if (local_v.size() != local_A.size()) {
    op.PreallocateWorkVectors();
    local_v = op.v;
    local_Av = op.Av;
  }

  Kokkos::parallel_for(
    "Operator_FaceCell::ApplyMatrixFreeOp",
    ncells_owned,
    KOKKOS_LAMBDA(const int& c){
      AmanziMesh::Entity_ID_View faces;
      mesh->cell_get_faces(c, faces);
      int nfaces = faces.size();

      WhetStone::DenseVector<Amanzi::DeviceOnlyMemorySpace> lv = getFromCSR<WhetStone::DenseVector>(local_v,c);
      WhetStone::DenseVector<Amanzi::DeviceOnlyMemorySpace> lAv = getFromCSR<WhetStone::DenseVector>(local_Av,c);
      WhetStone::DenseMatrix<Amanzi::DeviceOnlyMemorySpace> lA = getFromCSR<WhetStone::DenseMatrix>(local_A,c);

      for (int n = 0; n != nfaces; ++n) {
        lv(n) = Xf(0,faces[n]);
      }

      lA.Multiply(lv, lAv, false);

      for (int n = 0; n != nfaces; ++n) {
        Yf(0,faces[n]) += lAv(n);
      }
    }); 
  return 0;
}
#endif 

#if 0 
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
#endif 

#if 0 
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
#endif 

/* ******************************************************************
* Visit methods for symbolic assemble: FaceCell
****************************************************************** */
void Operator_FaceCell::SymbolicAssembleMatrixOp(const Op_Cell_FaceCell& op,
                                                 const SuperMap& map, GraphFE& graph,
                                                 int my_block_row, int my_block_col) const
{
  std::vector<int> lid_r(2*cell_max_faces + 1);
  std::vector<int> lid_c(2*cell_max_faces + 1);

  // ELEMENT: cell, DOFS: cell and face
  const auto face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const auto face_col_inds = map.GhostIndices(my_block_col, "face", 0);
  const auto cell_row_inds = map.GhostIndices(my_block_row, "cell", 0);
  const auto cell_col_inds = map.GhostIndices(my_block_col, "cell", 0);

  //const auto face_gh_map = map.ComponentGhostedMap(my_block_row, "face");
  //const auto cell_map = map.ComponentMap(my_block_row, "cell");
  
  //int ierr(0);
  AmanziMesh::Entity_ID_View faces;
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_faces(c, faces);
    int nfaces = faces.size();

    int k = 0;
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      lid_r[k] = face_row_inds[n];
      lid_c[k] = face_col_inds[n];
      ++k;
      //int face_dof_size = face_gh_map->ElementSize(f);
      //int first = face_gh_map->FirstPointInElement(f);

      //for (int m = 0; m != face_dof_size; ++m) {
      //  lid_r[k] = face_row_inds[first + m];
      //  lid_c[k] = face_col_inds[first + m];
      //  k++;
      //}
    }

    lid_r[k] = cell_row_inds[c]; 
    lid_c[k] = cell_col_inds[c];
    //int cell_dof_size = cell_map->ElementSize(c);
    //int first = cell_map->FirstPointInElement(c);
    //for (int m = 0; m != cell_dof_size; ++m) {
    //  lid_r[k] = cell_row_inds[first + m];
    //  lid_c[k] = cell_col_inds[first + m];
    //  k++;
    //}
    
    graph.insertLocalIndices(k, lid_r.data(), k, lid_c.data());
  }
  //AMANZI_ASSERT(!ierr);
}

#if 0 
/* ******************************************************************
* Visit methods for symbolic assemble: Face
****************************************************************** */
void Operator_FaceCell::SymbolicAssembleMatrixOp(const Op_Cell_Face& op,
                                                 const SuperMap& map, GraphFE& graph,
                                                 int my_block_row, int my_block_col) const
{
  std::vector<int> lid_r(2*cell_max_faces);
  std::vector<int> lid_c(2*cell_max_faces);

  // ELEMENT: cell, DOFS: face
  const auto face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const auto face_col_inds = map.GhostIndices(my_block_col, "face", 0);

  const auto face_gh_map = map.ComponentGhostedMap(my_block_row, "face");

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

      for (int m = 0; m != face_dof_size; ++m) {
        lid_r[k] = face_row_inds[first + m];
        lid_c[k] = face_col_inds[first + m];
        k++;
      }
    }
    ierr |= graph.insertLocalIndices(k, lid_r.data(), k, lid_c.data());
  }
  AMANZI_ASSERT(!ierr);
}
#endif 

#if 0 
/* ******************************************************************
* Visit methods for symbolic assemble: Surface
****************************************************************** */
void Operator_FaceCell::SymbolicAssembleMatrixOp(
    const Op_SurfaceCell_SurfaceCell& op,
    const SuperMap& map, GraphFE& graph,
    int my_block_row, int my_block_col) const
{
  int nsurf_cells = op.surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const std::vector<int>& face_col_inds = map.GhostIndices(my_block_col, "face", 0);

  int ierr = 0;
  for (int sc = 0; sc != nsurf_cells; ++sc) {
    int lid_r = face_row_inds[op.surf_mesh->entity_get_parent(AmanziMesh::CELL, sc)];
    int lid_c = face_col_inds[op.surf_mesh->entity_get_parent(AmanziMesh::CELL, sc)];
    ierr |= graph.insertLocalIndices(lid_r, 1, &lid_c);
  }
  AMANZI_ASSERT(!ierr);
}
#endif 

#if 0 
/* ******************************************************************
* Visit methods for symbolic assemble: Surface
****************************************************************** */
void Operator_FaceCell::SymbolicAssembleMatrixOp(
    const Op_SurfaceFace_SurfaceCell& op,
    const SuperMap& map, GraphFE& graph,
    int my_block_row, int my_block_col) const
{
  int nsurf_faces = op.surf_mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int lid_r[2];
  int lid_c[2];

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const std::vector<int>& face_col_inds = map.GhostIndices(my_block_col, "face", 0);

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

    ierr |= graph.insertLocalIndices(ncells, lid_r, ncells, lid_c);
  }
  AMANZI_ASSERT(!ierr);
}
#endif 

/* ******************************************************************
* Visit methods for assemble: FaceCell
****************************************************************** */
void Operator_FaceCell::AssembleMatrixOp(const Op_Cell_FaceCell& op,
                                         const SuperMap& map, MatrixFE& mat,
                                         int my_block_row, int my_block_col) const
{
  AMANZI_ASSERT(op.A.size() == ncells_owned);

  std::vector<int> lid_r(2*cell_max_faces + 1);
  std::vector<int> lid_c(2*cell_max_faces + 1);

  // ELEMENT: cell, DOFS: face and cell
  const auto face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const auto face_col_inds = map.GhostIndices(my_block_col, "face", 0);
  const auto cell_row_inds = map.GhostIndices(my_block_row, "cell", 0);
  const auto cell_col_inds = map.GhostIndices(my_block_col, "cell", 0);        

  const auto face_gh_map = map.ComponentGhostedMap(my_block_row, "face");
  const auto cell_map = map.ComponentMap(my_block_row, "cell");

  //int ierr(0);
  AmanziMesh::Entity_ID_View faces;
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_faces(c, faces);
    
    int nfaces = faces.size();
    int k = 0;
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      assert(false); 
      //int face_dof_size = face_gh_map->ElementSize(f);
      //int first = face_gh_map->FirstPointInElement(f);

      //for (int m = 0; m != face_dof_size; ++m) {
      //  lid_r[k] = face_row_inds[first + m];
      //  lid_c[k] = face_col_inds[first + m];
      //  k++;
      //}      
    }

    //int cell_dof_size = cell_map->ElementSize(c);
    //int first = cell_map->FirstPointInElement(c);
    //for (int m = 0; m != cell_dof_size; ++m) {
    //  lid_r[k] = cell_row_inds[first + m];
    //  lid_c[k] = cell_col_inds[first + m];
    //  k++;
    //}
    
    //ierr |= mat.sumIntoLocalValues(lid_r.data(), lid_c.data(), op.matrices[c]);
  }
  //AMANZI_ASSERT(!ierr);
}

#if 0 
/* ******************************************************************
* Visit methods for assemble: Face
****************************************************************** */
void Operator_FaceCell::AssembleMatrixOp(const Op_Cell_Face& op,
                                         const SuperMap& map, MatrixFE& mat,
                                         int my_block_row, int my_block_col) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

  std::vector<int> lid_r(cell_max_faces + 1);
  std::vector<int> lid_c(cell_max_faces + 1);

  // ELEMENT: cell, DOFS: face and cell
  const std::vector<int>& face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const std::vector<int>& face_col_inds = map.GhostIndices(my_block_col, "face", 0);

  Teuchos::RCP<const Epetra_BlockMap> face_gh_map = map.ComponentGhostedMap(my_block_row, "face");

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

      for (int m = 0; m != face_dof_size; ++m) {
        lid_r[k] = face_row_inds[first + m];
        lid_c[k] = face_col_inds[first + m];
        k++;
      }
    }    
    
    ierr |= mat.sumIntoLocalValues(lid_r.data(), lid_c.data(), op.matrices[c]);
  }
  AMANZI_ASSERT(!ierr);
}
#endif 

#if 0 
/* ******************************************************************
* Visit methods for assemble: Surface
****************************************************************** */
void Operator_FaceCell::AssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op,
                                         const SuperMap& map, MatrixFE& mat,
                                         int my_block_row, int my_block_col) const
{
  int nsurf_cells = op.surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const std::vector<int>& face_col_inds = map.GhostIndices(my_block_col, "face", 0);

  int ierr = 0;
  for (int sc = 0; sc != nsurf_cells; ++sc) {
    int lid_r = face_row_inds[op.surf_mesh->entity_get_parent(AmanziMesh::CELL, sc)];
    int lid_c = face_col_inds[op.surf_mesh->entity_get_parent(AmanziMesh::CELL, sc)];
    ierr |= mat.sumIntoLocalValues(lid_r, 1, &(*op.diag)[0][sc], &lid_c);
  }
  AMANZI_ASSERT(!ierr);
}
#endif 

#if 0 
void Operator_FaceCell::AssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op,
                                         const SuperMap& map, MatrixFE& mat,
                                         int my_block_row, int my_block_col) const
{
  int nsurf_faces = op.surf_mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  int lid_r[2];
  int lid_c[2];

  // ELEMENT: cell, DOFS: cell and face
  const std::vector<int>& face_row_inds = map.GhostIndices(my_block_row, "face", 0);
  const std::vector<int>& face_col_inds = map.GhostIndices(my_block_col, "face", 0);

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

    ierr |= mat.sumIntoLocalValues(lid_r, lid_c, op.matrices[sf]);
  }
  AMANZI_ASSERT(!ierr);
}
#endif 

}  // namespace Operators
}  // namespace Amanzi



