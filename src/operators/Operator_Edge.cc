/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)
*/

#include "DenseMatrix.hh"
#include "Op_Cell_Edge.hh"

#include "SuperMap.hh"
#include "GraphFE.hh"
#include "MatrixFE.hh"

#include "OperatorDefs.hh"
#include "Operator_Edge.hh"

/* ******************************************************************
Operator whose unknowns are Edge.

See Operator_Edge.hh for more detail.
****************************************************************** */

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Apply a source which may or may not have edge volume included already. 
****************************************************************** */
void Operator_Edge::UpdateRHS(const CompositeVector& source, bool volume_included)
{
  if (volume_included) {
    Operator::UpdateRHS(source);
  } else {
    rhs_->PutScalarGhosted(0.0);
    Epetra_MultiVector& rhs_e = *rhs_->ViewComponent("edge", true);
    const Epetra_MultiVector& source_e = *source.ViewComponent("edge", true);

    AmanziMesh::Entity_ID_List edges;

    for (int c = 0; c != ncells_owned; ++c) {
      mesh_->cell_get_edges(c, &edges);
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
int Operator_Edge::ApplyMatrixFreeOp(const Op_Cell_Edge& op,
                                     const CompositeVector& X, CompositeVector& Y) const
{
  ASSERT(op.matrices.size() == ncells_owned);

  X.ScatterMasterToGhosted();
  const Epetra_MultiVector& Xe = *X.ViewComponent("edge", true);
  Y.PutScalarGhosted(0.0);

  {
    Epetra_MultiVector& Ye = *Y.ViewComponent("edge", true);

    AmanziMesh::Entity_ID_List edges;
    for (int c = 0; c != ncells_owned; ++c) {
      mesh_->cell_get_edges(c, &edges);
      int nedges = edges.size();

      WhetStone::DenseVector v(nedges), av(nedges);
      for (int n = 0; n != nedges; ++n) {
        v(n) = Xe[0][edges[n]];
      }

      const WhetStone::DenseMatrix& Acell = op.matrices[c];
      Acell.Multiply(v, av, false);

      for (int n = 0; n != nedges; ++n) {
        Ye[0][edges[n]] += av(n);
      }
    } 
  }

  Y.GatherGhostedToMaster(Add);
  return 0;
}


/* ******************************************************************
* Visit methods for symbolic assemble.
* Apply the local matrices directly as schemas match.
****************************************************************** */
void Operator_Edge::SymbolicAssembleMatrixOp(const Op_Cell_Edge& op,
                                             const SuperMap& map, GraphFE& graph,
                                             int my_block_row, int my_block_col) const
{
  int lid_r[OPERATOR_MAX_EDGES];
  int lid_c[OPERATOR_MAX_EDGES];

  // ELEMENT: cell, DOFS: cell and edge
  const std::vector<int>& edge_row_inds = map.GhostIndices("edge", my_block_row);
  const std::vector<int>& edge_col_inds = map.GhostIndices("edge", my_block_col);

  int ierr(0);
  AmanziMesh::Entity_ID_List edges;
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_edges(c, &edges);
    int nedges = edges.size();

    for (int n = 0; n != nedges; ++n) {
      lid_r[n] = edge_row_inds[edges[n]];
      lid_c[n] = edge_col_inds[edges[n]];
    }
    ierr |= graph.InsertMyIndices(nedges, lid_r, nedges, lid_c);
  }
  ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble
* Apply the local matrices directly as schemas match.
****************************************************************** */
void Operator_Edge::AssembleMatrixOp(const Op_Cell_Edge& op,
                                     const SuperMap& map, MatrixFE& mat,
                                     int my_block_row, int my_block_col) const
{
  ASSERT(op.matrices.size() == ncells_owned);

  int lid_r[OPERATOR_MAX_EDGES];
  int lid_c[OPERATOR_MAX_EDGES];

  // ELEMENT: cell, DOFS: edge and cell
  const std::vector<int>& edge_row_inds = map.GhostIndices("edge", my_block_row);
  const std::vector<int>& edge_col_inds = map.GhostIndices("edge", my_block_col);

  int ierr(0);
  AmanziMesh::Entity_ID_List edges;
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_edges(c, &edges);
    int nedges = edges.size();

    for (int n = 0; n != nedges; ++n) {
      lid_r[n] = edge_row_inds[edges[n]];
      lid_c[n] = edge_col_inds[edges[n]];
    }

    ierr |= mat.SumIntoMyValues(lid_r, lid_c, op.matrices[c]);
  }
  ASSERT(!ierr);
}

}  // namespace Operators
}  // namespace Amanzi



