/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "AmanziTypes.hh"

#include "Op_MeshInjection.hh"
#include "SuperMap.hh"
#include "GraphFE.hh"
#include "MatrixFE.hh"
#include "Operator_Schema.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Op_MeshInjection represents one nonzero per row of diag's own index
* space (the "short" side). When transpose=false, that short side is
* mesh_row: diag_lid i is a row-local index, and injection->
* getGlobalElement(i) resolves (via mesh_col's CELL map) to the column.
* When transpose=true, the short side is mesh_col instead: diag_lid i is
* a COLUMN-local index, and injection->getGlobalElement(i) resolves (via
* mesh_row's CELL map) to the row.
****************************************************************** */
int
Operator_Schema::ApplyMatrixFreeOp(const Op_MeshInjection& op,
                                   const CompositeVector& X,
                                   CompositeVector& Y) const
{
  auto Xv = X.viewComponent("cell", false);
  auto Yv = Y.viewComponent("cell", false);
  const auto diag_v = op.diag->getLocalViewDevice(Tpetra::Access::ReadOnly);
  auto injection = op.injection;

  int n = op.diag->getLocalLength();
  if (!op.transpose) {
    // diag_lid == row_lid (the "short"/native side); col resolved via
    // injection + mesh_col's CELL map. Y[row_lid] += diag[row_lid] *
    // X[col_lid], matching AssembleMatrixOp's A[row_lid, col_lid] =
    // diag_v(row_lid, 0).
    const AmanziMesh::Mesh& mesh_col = op.get_col_mesh();
    const auto& col_cell_map = *mesh_col.getMap(AmanziMesh::Entity_kind::CELL, false);
    Kokkos::parallel_for(
      "Operator_Schema::ApplyMatrixFreeOp::MeshInjection", n, KOKKOS_LAMBDA(const int& row_lid) {
        auto col_gid = injection->getGlobalElement(row_lid);
        auto col_lid = col_cell_map.getLocalElement(col_gid);
        Kokkos::atomic_add(&Yv(row_lid, 0), diag_v(row_lid, 0) * Xv(col_lid, 0));
      });
  } else {
    // diag_lid == col_lid (the "short"/native side); row resolved via
    // injection + mesh_row's CELL map. Y[row_lid] += diag[col_lid] *
    // X[col_lid], matching AssembleMatrixOp's A[row_lid, col_lid] =
    // diag_v(col_lid, 0).
    const AmanziMesh::Mesh& mesh_row = op.get_row_mesh();
    const auto& row_cell_map = *mesh_row.getMap(AmanziMesh::Entity_kind::CELL, false);
    Kokkos::parallel_for(
      "Operator_Schema::ApplyMatrixFreeOp::MeshInjection", n, KOKKOS_LAMBDA(const int& col_lid) {
        auto row_gid = injection->getGlobalElement(col_lid);
        auto row_lid = row_cell_map.getLocalElement(row_gid);
        Kokkos::atomic_add(&Yv(row_lid, 0), diag_v(col_lid, 0) * Xv(col_lid, 0));
      });
  }
  return 0;
}


/* ******************************************************************
* Symbolic assembly: insert one (row,col) entry per diag-local index.
****************************************************************** */
void
Operator_Schema::SymbolicAssembleMatrixOp(const Op_MeshInjection& op,
                                          const SuperMap& map,
                                          GraphFE& graph,
                                          int my_block_row,
                                          int my_block_col) const
{
  const auto row_inds = map.viewGhostIndices<MemSpace_kind::HOST>(my_block_row, "cell", 0);
  const auto col_inds = map.viewGhostIndices<MemSpace_kind::HOST>(my_block_col, "cell", 0);
  auto injection = op.injection;
  int n = op.diag->getLocalLength();

  if (!op.transpose) {
    const AmanziMesh::Mesh& mesh_col = op.get_col_mesh();
    const auto& col_cell_map = *mesh_col.getMap(AmanziMesh::Entity_kind::CELL, false);
    for (int row_lid = 0; row_lid != n; ++row_lid) {
      auto col_gid = injection->getGlobalElement(row_lid);
      auto col_lid = col_cell_map.getLocalElement(col_gid);
      int row = row_inds[row_lid];
      int col = col_inds[col_lid];
      graph.insertLocalIndices(row, 1, &col);
    }
  } else {
    const AmanziMesh::Mesh& mesh_row = op.get_row_mesh();
    const auto& row_cell_map = *mesh_row.getMap(AmanziMesh::Entity_kind::CELL, false);
    for (int col_lid = 0; col_lid != n; ++col_lid) {
      auto row_gid = injection->getGlobalElement(col_lid);
      auto row_lid = row_cell_map.getLocalElement(row_gid);
      int row = row_inds[row_lid];
      int col = col_inds[col_lid];
      graph.insertLocalIndices(row, 1, &col);
    }
  }
}


/* ******************************************************************
* Assemble: insert one (row,col) scalar entry per diag-local index.
****************************************************************** */
void
Operator_Schema::AssembleMatrixOp(const Op_MeshInjection& op,
                                  const SuperMap& map,
                                  MatrixFE& mat,
                                  int my_block_row,
                                  int my_block_col) const
{
  const auto row_inds = map.viewGhostIndices<MemSpace_kind::HOST>(my_block_row, "cell", 0);
  const auto col_inds = map.viewGhostIndices<MemSpace_kind::HOST>(my_block_col, "cell", 0);
  const auto diag_v = op.diag->getLocalViewHost(Tpetra::Access::ReadOnly);
  auto injection = op.injection;
  int n = op.diag->getLocalLength();

  if (!op.transpose) {
    const AmanziMesh::Mesh& mesh_col = op.get_col_mesh();
    const auto& col_cell_map = *mesh_col.getMap(AmanziMesh::Entity_kind::CELL, false);
    for (int row_lid = 0; row_lid != n; ++row_lid) {
      auto col_gid = injection->getGlobalElement(row_lid);
      auto col_lid = col_cell_map.getLocalElement(col_gid);
      int row = row_inds[row_lid];
      int col = col_inds[col_lid];
      double val = diag_v(row_lid, 0);
      mat.sumIntoLocalValues(row, 1, &val, &col);
    }
  } else {
    const AmanziMesh::Mesh& mesh_row = op.get_row_mesh();
    const auto& row_cell_map = *mesh_row.getMap(AmanziMesh::Entity_kind::CELL, false);
    for (int col_lid = 0; col_lid != n; ++col_lid) {
      auto row_gid = injection->getGlobalElement(col_lid);
      auto row_lid = row_cell_map.getLocalElement(row_gid);
      int row = row_inds[row_lid];
      int col = col_inds[col_lid];
      double val = diag_v(col_lid, 0);
      mat.sumIntoLocalValues(row, 1, &val, &col);
    }
  }
}

} // namespace Operators
} // namespace Amanzi
