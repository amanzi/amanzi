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

  Operator whose unknowns are Edge.
*/

#include "DenseMatrix.hh"
#include "Op_Edge_Edge.hh"
#include "Op_Cell_Edge.hh"

#include "SuperMap.hh"
#include "GraphFE.hh"
#include "MatrixFE.hh"

#include "OperatorDefs.hh"
#include "Operator_Edge.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Apply a source which may or may not have edge volume included already.
****************************************************************** */
void
Operator_Edge::UpdateRHS(const CompositeVector& source, bool volume_included)
{
  if (volume_included) {
    Operator::UpdateRHS(source);
  } else {
    rhs_->PutScalarGhosted(0.0);
    Epetra_MultiVector& rhs_e = *rhs_->ViewComponent("edge", true);
    const Epetra_MultiVector& source_e = *source.ViewComponent("edge", true);

    for (int c = 0; c != ncells_owned; ++c) {
      const auto& edges = mesh_->cell_get_edges(c);
      int nedges = edges.size();

      double volume = mesh_->cell_volume(c);
      for (int n = 0; n < nedges; ++n) {
        int e = edges[n];
        rhs_e[0][e] += source_e[0][e] * volume / nedges;
      }
    }
    rhs_->GatherGhostedToMaster("edge");
  }
}


/* ******************************************************************
* Apply the local matrices directly as schemas match.
****************************************************************** */
int
Operator_Edge::ApplyMatrixFreeOp(const Op_Cell_Edge& op,
                                 const CompositeVector& X,
                                 CompositeVector& Y) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);
  const Epetra_MultiVector& Xe = *X.ViewComponent("edge", true);
  {
    Epetra_MultiVector& Ye = *Y.ViewComponent("edge", true);

    AmanziMesh::Entity_ID_List edges;
    for (int c = 0; c != ncells_owned; ++c) {
      mesh_->cell_get_edges(c, &edges);
      int nedges = edges.size();

      WhetStone::DenseVector v(nedges), av(nedges);
      for (int n = 0; n != nedges; ++n) { v(n) = Xe[0][edges[n]]; }

      const WhetStone::DenseMatrix& Acell = op.matrices[c];
      Acell.Multiply(v, av, false);

      for (int n = 0; n != nedges; ++n) { Ye[0][edges[n]] += av(n); }
    }
  }
  return 0;
}


/* ******************************************************************
* Visit methods for Apply.
* Apply the local matrices directly as schema is a subset of
* assembled schema.
****************************************************************** */
int
Operator_Edge::ApplyMatrixFreeOp(const Op_Edge_Edge& op,
                                 const CompositeVector& X,
                                 CompositeVector& Y) const
{
  const Epetra_MultiVector& Xc = *X.ViewComponent("edge");
  Epetra_MultiVector& Yc = *Y.ViewComponent("edge");

  for (int k = 0; k != Xc.NumVectors(); ++k) {
    for (int e = 0; e != nedges_owned; ++e) { Yc[k][e] += Xc[k][e] * (*op.diag)[k][e]; }
  }
  return 0;
}


/* ******************************************************************
* Visit methods for symbolic assemble.
* Apply the local matrices directly as schemas match.
****************************************************************** */
void
Operator_Edge::SymbolicAssembleMatrixOp(const Op_Cell_Edge& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const
{
  std::vector<int> lid_r(cell_max_edges);
  std::vector<int> lid_c(cell_max_edges);

  // ELEMENT: cell, DOFS: cell and edge
  const std::vector<int>& edge_row_inds = map.GhostIndices(my_block_row, "edge", 0);
  const std::vector<int>& edge_col_inds = map.GhostIndices(my_block_col, "edge", 0);

  int ierr(0);
  AmanziMesh::Entity_ID_List edges;
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_edges(c, &edges);
    int nedges = edges.size();

    for (int n = 0; n != nedges; ++n) {
      lid_r[n] = edge_row_inds[edges[n]];
      lid_c[n] = edge_col_inds[edges[n]];
    }
    ierr |= graph.InsertMyIndices(nedges, lid_r.data(), nedges, lid_c.data());
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for symbolic assemble.
* Insert the diagonal on edges
****************************************************************** */
void
Operator_Edge::SymbolicAssembleMatrixOp(const Op_Edge_Edge& op,
                                        const SuperMap& map,
                                        GraphFE& graph,
                                        int my_block_row,
                                        int my_block_col) const
{
  const std::vector<int>& edge_row_inds = map.GhostIndices(my_block_row, "edge", 0);
  const std::vector<int>& edge_col_inds = map.GhostIndices(my_block_col, "edge", 0);

  int ierr(0);
  for (int e = 0; e != nedges_owned; ++e) {
    int row = edge_row_inds[e];
    int col = edge_col_inds[e];

    ierr |= graph.InsertMyIndices(row, 1, &col);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble
* Apply the local matrices directly as schemas match.
****************************************************************** */
void
Operator_Edge::AssembleMatrixOp(const Op_Cell_Edge& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

  std::vector<int> lid_r(cell_max_edges);
  std::vector<int> lid_c(cell_max_edges);

  // ELEMENT: cell, DOFS: edge and cell
  const std::vector<int>& edge_row_inds = map.GhostIndices(my_block_row, "edge", 0);
  const std::vector<int>& edge_col_inds = map.GhostIndices(my_block_col, "edge", 0);

  int ierr(0);
  AmanziMesh::Entity_ID_List edges;
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_edges(c, &edges);
    int nedges = edges.size();

    for (int n = 0; n != nedges; ++n) {
      lid_r[n] = edge_row_inds[edges[n]];
      lid_c[n] = edge_col_inds[edges[n]];
    }

    ierr |= mat.SumIntoMyValues(lid_r.data(), lid_c.data(), op.matrices[c]);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble
* Insert each diagonal values for edges.
****************************************************************** */
void
Operator_Edge::AssembleMatrixOp(const Op_Edge_Edge& op,
                                const SuperMap& map,
                                MatrixFE& mat,
                                int my_block_row,
                                int my_block_col) const
{
  AMANZI_ASSERT(op.diag->NumVectors() == 1);

  const std::vector<int>& edge_row_inds = map.GhostIndices(my_block_row, "edge", 0);
  const std::vector<int>& edge_col_inds = map.GhostIndices(my_block_col, "edge", 0);

  int ierr(0);
  for (int e = 0; e != nedges_owned; ++e) {
    int row = edge_row_inds[e];
    int col = edge_col_inds[e];

    ierr |= mat.SumIntoMyValues(row, 1, &(*op.diag)[0][e], &col);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Operator>
Operator_Edge::Clone() const
{
  return Teuchos::rcp(new Operator_Edge(*this));
}

} // namespace Operators
} // namespace Amanzi
