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
#include "Op_Cell_Cell.hh"
#include "Op_Face_Cell.hh"
#include "OperatorUtils.hh"

#include "SuperMap.hh"
#include "GraphFE.hh"
#include "MatrixFE.hh"
#include "Operator_Cell.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
 * Apply a source which may or may not have cell volume included already.
 ****************************************************************** */
void
Operator_Cell::UpdateRHS(const CompositeVector& source, bool volume_included)
{
  if (volume_included) {
    Operator::UpdateRHS(source);
  } else {
    auto rhs_v = rhs_->ViewComponent("cell", false);
    auto source_v = source.ViewComponent("cell", false);
    const AmanziMesh::Mesh* mesh_p = mesh_.get();
    Kokkos::parallel_for(
        ncells_owned,
        KOKKOS_LAMBDA(const int c) {
          rhs_v(c,0) += source_v(c,0) * mesh_p->cell_volume(c);
        });
  }
}


/* ******************************************************************
 * Visit methods for Apply.
 * Apply the local matrices directly as schema is a subset of
 * assembled schema.
 ****************************************************************** */
int
Operator_Cell::ApplyMatrixFreeOp(const Op_Cell_Cell& op,
                                 const CompositeVector& X,
                                 CompositeVector& Y) const
{
  auto X_v = X.ViewComponent("cell");
  auto Y_v = Y.ViewComponent("cell");
  auto& op_v = op.data;

  // currently hard-coded in Op_Cell_Cell, though this could/should be relaxed.
  AMANZI_ASSERT(X_v.extent(1) == 1);
  AMANZI_ASSERT(Y_v.extent(1) == 1);
  
  Kokkos::MDRangePolicy<Kokkos::Rank<2>> policy({0,0},
          {X_v.extent(0), X_v.extent(1)});
  Kokkos::parallel_for(policy,
                       KOKKOS_LAMBDA(const int i, const int j) {
                         Y_v(i,j) += X_v(i,j) * op_v(i,j);
                       });
  return 0;
}


/* ******************************************************************
 * Apply the local matrices directly as schema is a subset of
 * assembled schema
 ****************************************************************** */
int
Operator_Cell::ApplyMatrixFreeOp(const Op_Face_Cell& op,
                                 const CompositeVector& X,
                                 CompositeVector& Y) const
{
  AMANZI_ASSERT(op.data.extent(0) == nfaces_owned);

  X.ScatterMasterToGhosted();
  Y.putScalarGhosted(0.);

  {
    auto X_v = X.ViewComponent("cell", true);
    auto Y_v = Y.ViewComponent("cell", true);
    auto mesh_p = mesh_.get();
    Kokkos::View<const double**> op_v = op.data;

    Kokkos::parallel_for(
        "local_mat_mult_element_loop",
        op_v.extent(0),
        KOKKOS_LAMBDA(const int f) {
          AmanziMesh::Entity_ID_View cells;
          mesh_p->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
          auto num_cols = cells.extent(0);
          auto num_rows = cells.extent(0);
          for(int i = 0 ; i < num_rows; ++i){
            double result = 0.0; 
            for(int j = 0 ; j < cells.extent(0); ++j){
              result += op_v(f,j+i*num_cols)*X_v(cells(j),0); 
            }
            Kokkos::atomic_add(&Y_v(cells(i), 0), result);
          }
        });
  }

  Y.GatherGhostedToMaster("cell", Tpetra::ADD);
  return 0;
}


// /* ******************************************************************
//  * Visit methods for symbolic assemble.
//  * Insert the diagonal on cells
//  ****************************************************************** */
// void
// Operator_Cell::SymbolicAssembleMatrixOp(const Op_Cell_Cell& op,
//                                         const SuperMap& map, GraphFE& graph,
//                                         int my_block_row,
//                                         int my_block_col) const
// {
//   const std::vector<int>& cell_row_inds =
//     map.GhostIndices(my_block_row, "cell", 0);
//   const std::vector<int>& cell_col_inds =
//     map.GhostIndices(my_block_col, "cell", 0);

//   int ierr(0);
//   for (int c = 0; c != ncells_owned; ++c) {
//     int row = cell_row_inds[c];
//     int col = cell_col_inds[c];

//     ierr |= graph.InsertMyIndices(row, 1, &col);
//   }
//   AMANZI_ASSERT(!ierr);
// }


// /* ******************************************************************
//  * Insert each cells neighboring cells.
//  ****************************************************************** */
// void
// Operator_Cell::SymbolicAssembleMatrixOp(const Op_Face_Cell& op,
//                                         const SuperMap& map, GraphFE& graph,
//                                         int my_block_row,
//                                         int my_block_col) const
// {
//   // ELEMENT: face, DOF: cell
//   int lid_r[2];
//   int lid_c[2];

//   const std::vector<int>& cell_row_inds =
//     map.GhostIndices(my_block_row, "cell", 0);
//   const std::vector<int>& cell_col_inds =
//     map.GhostIndices(my_block_col, "cell", 0);

//   int ierr(0);
//   AmanziMesh::Entity_ID_List cells;
//   for (int f = 0; f != nfaces_owned; ++f) {
//     mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

//     int ncells = cells.size();
//     for (int n = 0; n != ncells; ++n) {
//       lid_r[n] = cell_row_inds[cells[n]];
//       lid_c[n] = cell_col_inds[cells[n]];
//     }

//     ierr |= graph.InsertMyIndices(ncells, lid_r, ncells, lid_c);
//   }

//   AMANZI_ASSERT(!ierr);
// }


// /* ******************************************************************
//  * Visit methods for assemble
//  * Insert each cells neighboring cells.
//  ****************************************************************** */
// void
// Operator_Cell::AssembleMatrixOp(const Op_Cell_Cell& op, const SuperMap& map,
//                                 MatrixFE& mat, int my_block_row,
//                                 int my_block_col) const
// {
//   AMANZI_ASSERT(op.diag->getNumVectors() == 1);

//   const std::vector<int>& cell_row_inds =
//     map.GhostIndices(my_block_row, "cell", 0);
//   const std::vector<int>& cell_col_inds =
//     map.GhostIndices(my_block_col, "cell", 0);

//   int ierr(0);
//   for (int c = 0; c != ncells_owned; ++c) {
//     int row = cell_row_inds[c];
//     int col = cell_col_inds[c];

//     ierr |= mat.SumIntoMyValues(row, 1, &(*op.diag)[0][c], &col);
//   }
//   AMANZI_ASSERT(!ierr);
// }


// void
// Operator_Cell::AssembleMatrixOp(const Op_Face_Cell& op, const SuperMap& map,
//                                 MatrixFE& mat, int my_block_row,
//                                 int my_block_col) const
// {
//   AMANZI_ASSERT(op.matrices.size() == nfaces_owned);

//   // ELEMENT: face, DOF: cell
//   int lid_r[2];
//   int lid_c[2];

//   const std::vector<int>& cell_row_inds =
//     map.GhostIndices(my_block_row, "cell", 0);
//   const std::vector<int>& cell_col_inds =
//     map.GhostIndices(my_block_col, "cell", 0);

//   int ierr(0);
//   AmanziMesh::Entity_ID_List cells;
//   for (int f = 0; f != nfaces_owned; ++f) {
//     mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

//     int ncells = cells.size();
//     for (int n = 0; n != ncells; ++n) {
//       lid_r[n] = cell_row_inds[cells[n]];
//       lid_c[n] = cell_col_inds[cells[n]];
//     }

//     ierr |= mat.SumIntoMyValues(lid_r, lid_c, op.matrices[f]);
//     AMANZI_ASSERT(!ierr);
//   }
//   AMANZI_ASSERT(!ierr);
// }

} // namespace Operators
} // namespace Amanzi
