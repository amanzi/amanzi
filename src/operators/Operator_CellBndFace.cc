/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Daniil Svyatsky(dasvyat@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "DenseMatrix.hh"
#include "Op_Cell_Cell.hh"
#include "Op_Face_CellBndFace.hh"
#include "Op_SurfaceCell_SurfaceCell.hh"
#include "Op_SurfaceFace_SurfaceCell.hh"


#include "SuperMap.hh"
#include "GraphFE.hh"
#include "MatrixFE.hh"
#include "Operator_CellBndFace.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
 * Apply the local matrices directly as schema is a subset of
 * assembled schema
 ****************************************************************** */
int
Operator_CellBndFace::ApplyMatrixFreeOp(const Op_Face_CellBndFace& op,
                                        const CompositeVector& X,
                                        CompositeVector& Y) const
{
  AMANZI_ASSERT(op.data.extent(0) == nfaces_owned);

  X.ScatterMasterToGhosted();
  Y.putScalarGhosted(0.);

  // Stage views for kernel
  {
    // -- views into X and Y
    auto X_c = X.ViewComponent("cell", true);
    auto X_bf = X.ViewComponent("boundary_face", true);
    auto Y_c = Y.ViewComponent("cell", true);
    auto Y_bf = Y.ViewComponent("boundary_face", true);

    // -- views into our matrix and mesh
    auto mesh_p = mesh_.get();
    Kokkos::View<const double**> op_v = op.data;

    Kokkos::parallel_for(
        "local_mat_mult_element_loop",
        op_v.extent(0),
        KOKKOS_LAMBDA(const int f) {
          AmanziMesh::Entity_ID_View cells;
          mesh_p->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
          int ncells = cells.size();

          Kokkos::atomic_add(&Y_c(cells(0),0), op_v(f,0) * X_c(cells(0),0));

          if (ncells == 2) {
            Kokkos::atomic_add(&Y_c(cells(0),0), op_v(f,1) * X_c(cells(1),0));
            Kokkos::atomic_add(&Y_c(cells(1),0), op_v(f,2) * X_c(cells(0),0));
            Kokkos::atomic_add(&Y_c(cells(1),0), op_v(f,3) * X_c(cells(1),0));
          }
        });

    // now do boundary faces
    // -- views from bf to f
    // The importer gets us the face LID of a given boundary face LID
    // through two parts -- the first import_same_face are bf = f, then the
    // remainder are permutations.  Note this is guaranteed because we
    // know, by definition, that the faces and boundary faces are owned by
    // the same process.        
    auto import_same_face = mesh_->exterior_face_importer()->getNumSameIDs();
    auto import_permute_face =
        mesh_->exterior_face_importer()
        ->getPermuteFromLIDs_dv().view<AmanziDefaultDevice>();

    Kokkos::parallel_for(
        "local_mat_mult_element_loop_bf_same",
        import_same_face,
        KOKKOS_LAMBDA(const int& bf) {
          AmanziMesh::Entity_ID_View cells;
          mesh_p->face_get_cells(bf, AmanziMesh::Parallel_type::ALL, cells);

          Kokkos::atomic_add(&Y_c(cells(0),0), op_v(bf,1) * X_bf(bf,0));
          Kokkos::atomic_add(&Y_bf(bf,0), op_v(bf,2) * X_c(cells(0),0));
          Kokkos::atomic_add(&Y_bf(bf,0), op_v(bf,3) * X_bf(bf,0));
        });

    Kokkos::parallel_for(
        "local_mat_mult_element_loop_bf_permute",
        import_permute_face.extent(0),
        KOKKOS_LAMBDA(const int& bf_offset) {
          int bf = import_same_face + bf_offset;
          int f = import_permute_face(bf_offset);
          AmanziMesh::Entity_ID_View cells;
          mesh_p->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);

          Kokkos::atomic_add(&Y_c(cells(0),0), op_v(f,1) * X_bf(bf,0));
          Kokkos::atomic_add(&Y_bf(bf,0), op_v(f,2) * X_c(cells(0),0));
          Kokkos::atomic_add(&Y_bf(bf,0), op_v(f,3) * X_bf(bf,0));
        });
  }
  
  Y.GatherGhostedToMaster("cell", Tpetra::ADD);
  Y.GatherGhostedToMaster("boundary_face", Tpetra::ADD);
  return 0;
}


// /* ******************************************************************
//  * Insert each cells neighboring cells.
//  ****************************************************************** */
// void
// Operator_CellBndFace::SymbolicAssembleMatrixOp(const Op_Face_CellBndFace& op,
//                                                const SuperMap& map,
//                                                GraphFE& graph, int my_block_row,
//                                                int my_block_col) const
// {
//   // ELEMENT: face, DOF: cell, bnd_face
//   int lid_r[2];
//   int lid_c[2];
//   const std::vector<int>& cell_row_inds =
//     map.GhostIndices(my_block_row, "cell", 0);
//   const std::vector<int>& cell_col_inds =
//     map.GhostIndices(my_block_col, "cell", 0);
//   const std::vector<int>& bndface_row_inds =
//     map.GhostIndices(my_block_row, "boundary_face", 0);
//   const std::vector<int>& bndface_col_inds =
//     map.GhostIndices(my_block_col, "boundary_face", 0);

//   int ierr(0);
//   AmanziMesh::Entity_ID_List cells;
//   for (int f = 0; f != nfaces_owned; ++f) {
//     mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

//     int ncells = cells.size();
//     if (ncells == 2) {
//       for (int n = 0; n != ncells; ++n) {
//         lid_r[n] = cell_row_inds[cells[n]];
//         lid_c[n] = cell_col_inds[cells[n]];
//       }
//     } else if (ncells == 1) {
//       lid_r[0] = cell_row_inds[cells[0]];
//       lid_c[0] = cell_col_inds[cells[0]];
//       int bf = mesh_->exterior_face_map(false).getLocalElement(
//         mesh_->face_map(false).getGlobalElement(f));
//       lid_r[1] = bndface_row_inds[bf];
//       lid_c[1] = bndface_col_inds[bf];
//     }


//     ierr |= graph.InsertMyIndices(2, lid_r, 2, lid_c);
//   }
//   AMANZI_ASSERT(!ierr);
// }


// void
// Operator_CellBndFace::AssembleMatrixOp(const Op_Face_CellBndFace& op,
//                                        const SuperMap& map, MatrixFE& mat,
//                                        int my_block_row, int my_block_col) const
// {
//   AMANZI_ASSERT(op.matrices.size() == nfaces_owned);

//   // ELEMENT: face, DOF: cell,  bnd_face
//   int lid_r[2];
//   int lid_c[2];
//   const std::vector<int>& cell_row_inds =
//     map.GhostIndices(my_block_row, "cell", 0);
//   const std::vector<int>& cell_col_inds =
//     map.GhostIndices(my_block_col, "cell", 0);
//   const std::vector<int>& bndface_row_inds =
//     map.GhostIndices(my_block_row, "boundary_face", 0);
//   const std::vector<int>& bndface_col_inds =
//     map.GhostIndices(my_block_col, "boundary_face", 0);

//   int ierr(0);
//   AmanziMesh::Entity_ID_List cells;
//   for (int f = 0; f != nfaces_owned; ++f) {
//     mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

//     int ncells = cells.size();
//     if (ncells == 2) {
//       for (int n = 0; n != ncells; ++n) {
//         lid_r[n] = cell_row_inds[cells[n]];
//         lid_c[n] = cell_col_inds[cells[n]];
//       }
//     } else if (ncells == 1) {
//       lid_r[0] = cell_row_inds[cells[0]];
//       lid_c[0] = cell_col_inds[cells[0]];
//       int bf = mesh_->exterior_face_map(false).getLocalElement(
//         mesh_->face_map(false).getGlobalElement(f));
//       lid_r[1] = bndface_row_inds[bf];
//       lid_c[1] = bndface_col_inds[bf];
//     }

//     ierr |= mat.SumIntoMyValues(lid_r, lid_c, op.matrices[f]);
//     AMANZI_ASSERT(!ierr);
//   }
//   AMANZI_ASSERT(!ierr);
// }

// int
// Operator_CellBndFace::ApplyMatrixFreeOp(const Op_SurfaceCell_SurfaceCell& op,
//                                         const CompositeVector& X,
//                                         CompositeVector& Y) const
// {
//   int nsurf_cells = op.surf_mesh->num_entities(
//     AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
//   AMANZI_ASSERT(op.diag->getLocalLength() == nsurf_cells);

//   const Epetra_MultiVector& Xf = *X.ViewComponent("boundary_face", false);
//   Epetra_MultiVector& Yf = *Y.ViewComponent("boundary_face", false);
//   for (int sc = 0; sc != nsurf_cells; ++sc) {
//     int f = op.surf_mesh->entity_get_parent(AmanziMesh::CELL, sc);
//     int bf = mesh_->exterior_face_map(false).getLocalElement(
//       mesh_->face_map(false).getGlobalElement(f));
//     Yf[0][bf] += (*op.diag)[0][sc] * Xf[0][bf];
//   }
//   return 0;
// }

// /* ******************************************************************
//  * visit method for apply surface cells into subsurface faces
//  ****************************************************************** */
// int
// Operator_CellBndFace::ApplyMatrixFreeOp(const Op_SurfaceFace_SurfaceCell& op,
//                                         const CompositeVector& X,
//                                         CompositeVector& Y) const
// {
//   int nsurf_faces = op.surf_mesh->num_entities(
//     AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
//   AMANZI_ASSERT(op.matrices.size() == nsurf_faces);

//   X.ScatterMasterToGhosted();
//   const Epetra_MultiVector& Xf = *X.ViewComponent("boundary_face", true);

//   Y.putScalarGhosted(0.);
//   Epetra_MultiVector& Yf = *Y.ViewComponent("boundary_face", true);

//   AmanziMesh::Entity_ID_List cells;
//   for (int sf = 0; sf != nsurf_faces; ++sf) {
//     op.surf_mesh->face_get_cells(sf, AmanziMesh::Parallel_type::ALL, &cells);
//     int ncells = cells.size();

//     WhetStone::DenseVector v(ncells), av(ncells);
//     for (int n = 0; n != ncells; ++n) {
//       int f = op.surf_mesh->entity_get_parent(AmanziMesh::CELL, cells[n]);
//       int bf = mesh_->exterior_face_map(true).getLocalElement(
//         mesh_->face_map(true).getGlobalElement(f));
//       v(n) = Xf[0][bf];
//     }

//     const WhetStone::DenseMatrix& Aface = op.matrices[sf];
//     Aface.elementWiseMultiply(v, av, false);

//     for (int n = 0; n != ncells; ++n) {
//       int f = op.surf_mesh->entity_get_parent(AmanziMesh::CELL, cells[n]);
//       int bf = mesh_->exterior_face_map(true).getLocalElement(
//         mesh_->face_map(true).getGlobalElement(f));
//       Yf[0][bf] += av(n);
//     }
//   }
//   Y.GatherGhostedToMaster("boundary_face");
//   return 0;
// }


// void
// Operator_CellBndFace::SymbolicAssembleMatrixOp(
//   const Op_SurfaceCell_SurfaceCell& op, const SuperMap& map, GraphFE& graph,
//   int my_block_row, int my_block_col) const
// {
//   int nsurf_cells = op.surf_mesh->num_entities(
//     AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

//   // ELEMENT: cell, DOFS: cell and face
//   const std::vector<int>& face_row_inds =
//     map.GhostIndices(my_block_row, "boundary_face", 0);
//   const std::vector<int>& face_col_inds =
//     map.GhostIndices(my_block_col, "boundary_face", 0);

//   int ierr = 0;
//   for (int sc = 0; sc != nsurf_cells; ++sc) {
//     int f = op.surf_mesh->entity_get_parent(AmanziMesh::CELL, sc);
//     int bf = mesh_->exterior_face_map(true).getLocalElement(
//       mesh_->face_map(false).getGlobalElement(f));
//     int lid_r = face_row_inds[bf];
//     int lid_c = face_col_inds[bf];
//     ierr |= graph.InsertMyIndices(lid_r, 1, &lid_c);
//   }
//   AMANZI_ASSERT(!ierr);
// }


// void
// Operator_CellBndFace::SymbolicAssembleMatrixOp(
//   const Op_SurfaceFace_SurfaceCell& op, const SuperMap& map, GraphFE& graph,
//   int my_block_row, int my_block_col) const
// {
//   int nsurf_faces = op.surf_mesh->num_entities(
//     AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
//   int lid_r[2];
//   int lid_c[2];

//   // ELEMENT: cell, DOFS: cell and face
//   const std::vector<int>& face_row_inds =
//     map.GhostIndices(my_block_row, "boundary_face", 0);
//   const std::vector<int>& face_col_inds =
//     map.GhostIndices(my_block_col, "boundary_face", 0);

//   int ierr = 0;
//   AmanziMesh::Entity_ID_List cells;
//   for (int sf = 0; sf != nsurf_faces; ++sf) {
//     op.surf_mesh->face_get_cells(sf, AmanziMesh::Parallel_type::ALL, &cells);
//     int ncells = cells.size();
//     for (int n = 0; n != ncells; ++n) {
//       int f = op.surf_mesh->entity_get_parent(AmanziMesh::CELL, cells[n]);
//       int bf = mesh_->exterior_face_map(true).getLocalElement(
//         mesh_->face_map(true).getGlobalElement(f));
//       lid_r[n] = face_row_inds[bf];
//       lid_c[n] = face_col_inds[bf];
//     }

//     ierr |= graph.InsertMyIndices(ncells, lid_r, ncells, lid_c);
//   }
//   AMANZI_ASSERT(!ierr);
//   //   exit(0);
// }


// void
// Operator_CellBndFace::AssembleMatrixOp(const Op_SurfaceCell_SurfaceCell& op,
//                                        const SuperMap& map, MatrixFE& mat,
//                                        int my_block_row, int my_block_col) const
// {
//   int nsurf_cells = op.surf_mesh->num_entities(
//     AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

//   // ELEMENT: cell, DOFS: cell and face
//   const std::vector<int>& face_row_inds =
//     map.GhostIndices(my_block_row, "boundary_face", 0);
//   const std::vector<int>& face_col_inds =
//     map.GhostIndices(my_block_col, "boundary_face", 0);

//   int ierr = 0;
//   for (int sc = 0; sc != nsurf_cells; ++sc) {
//     int f = op.surf_mesh->entity_get_parent(AmanziMesh::CELL, sc);
//     int bf = mesh_->exterior_face_map(true).getLocalElement(
//       mesh_->face_map(true).getGlobalElement(f));

//     int lid_r = face_row_inds[bf];
//     int lid_c = face_col_inds[bf];
//     ierr |= mat.SumIntoMyValues(lid_r, 1, &(*op.diag)[0][sc], &lid_c);
//   }
//   AMANZI_ASSERT(!ierr);
// }


// void
// Operator_CellBndFace::AssembleMatrixOp(const Op_SurfaceFace_SurfaceCell& op,
//                                        const SuperMap& map, MatrixFE& mat,
//                                        int my_block_row, int my_block_col) const
// {
//   int nsurf_faces = op.surf_mesh->num_entities(
//     AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
//   int lid_r[2];
//   int lid_c[2];

//   // ELEMENT: cell, DOFS: cell and face
//   const std::vector<int>& face_row_inds =
//     map.GhostIndices(my_block_row, "boundary_face", 0);
//   const std::vector<int>& face_col_inds =
//     map.GhostIndices(my_block_col, "boundary_face", 0);

//   int ierr = 0;
//   AmanziMesh::Entity_ID_List cells;
//   for (int sf = 0; sf != nsurf_faces; ++sf) {
//     op.surf_mesh->face_get_cells(sf, AmanziMesh::Parallel_type::ALL, &cells);
//     int ncells = cells.size();
//     for (int n = 0; n != ncells; ++n) {
//       int f = op.surf_mesh->entity_get_parent(AmanziMesh::CELL, cells[n]);
//       int bf = mesh_->exterior_face_map(true).getLocalElement(
//         mesh_->face_map(true).getGlobalElement(f));
//       lid_r[n] = face_row_inds[bf];
//       lid_c[n] = face_col_inds[bf];
//     }

//     ierr |= mat.SumIntoMyValues(lid_r, lid_c, op.matrices[sf]);
//   }

//   AMANZI_ASSERT(!ierr);
// }


} // namespace Operators
} // namespace Amanzi
