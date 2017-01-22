/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  Operator whose unknowns are NODEs.
*/

#include "DenseMatrix.hh"
#include "GraphFE.hh"
#include "MatrixFE.hh"
#include "SuperMap.hh"

#include "Op_Cell_Node.hh"
#include "Op_Node_Node.hh"
#include "OperatorDefs.hh"
#include "Operator_Node.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Apply a source which may or may not have node volume included already. 
****************************************************************** */
void Operator_Node::UpdateRHS(const CompositeVector& source, bool volume_included)
{
  if (volume_included) {
    Operator::UpdateRHS(source);
  } else {
    rhs_->PutScalarGhosted(0.0);
    Epetra_MultiVector& rhs_v = *rhs_->ViewComponent("node", true);
    const Epetra_MultiVector& source_v = *source.ViewComponent("node", true);

    AmanziMesh::Entity_ID_List nodes;

    for (int c = 0; c != ncells_owned; ++c) {
      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      double volume = mesh_->cell_volume(c);
      for (int n = 0; n < nnodes; ++n) {
        int v = nodes[n];
        rhs_v[0][v] += source_v[0][v] * volume / nnodes;
      }
    }
    rhs_->GatherGhostedToMaster("node");
  }
}


/* ******************************************************************
* Apply the local matrices directly as schemas match.
****************************************************************** */
int Operator_Node::ApplyMatrixFreeOp(const Op_Cell_Node& op,
                                     const CompositeVector& X, CompositeVector& Y) const
{
  ASSERT(op.matrices.size() == ncells_owned);

  X.ScatterMasterToGhosted();
  const Epetra_MultiVector& Xn = *X.ViewComponent("node", true);
  Y.PutScalarGhosted(0.0);

  {
    Epetra_MultiVector& Yn = *Y.ViewComponent("node", true);

    AmanziMesh::Entity_ID_List nodes;
    for (int c = 0; c != ncells_owned; ++c) {
      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      WhetStone::DenseVector v(nnodes), av(nnodes);
      for (int n = 0; n != nnodes; ++n) {
        v(n) = Xn[0][nodes[n]];
      }

      const WhetStone::DenseMatrix& Acell = op.matrices[c];
      Acell.Multiply(v, av, false);

      for (int n = 0; n != nnodes; ++n) {
        Yn[0][nodes[n]] += av(n);
      }
    } 
  }

  Y.GatherGhostedToMaster(Add);
  return 0;
}


/* ******************************************************************
* Apply the local matrices directly as schemas match.
****************************************************************** */
int Operator_Node::ApplyMatrixFreeOp(const Op_Node_Node& op,
                                     const CompositeVector& X, CompositeVector& Y) const
{
  ASSERT(op.vals.size() == nnodes_owned);
  const Epetra_MultiVector& Xn = *X.ViewComponent("node");
  Epetra_MultiVector& Yn = *Y.ViewComponent("node");

  for (int i = 0; i < Xn.NumVectors(); ++i) {
    for (int v = 0; v != nnodes_owned; ++v) {
      Yn[i][v] += Xn[i][v] * op.vals[v];
    }
  }
  return 0;
}

/* ******************************************************************
* Visit methods for symbolic assemble.
* Apply the local matrices directly as schemas match.
****************************************************************** */
void Operator_Node::SymbolicAssembleMatrixOp(const Op_Cell_Node& op,
                                             const SuperMap& map, GraphFE& graph,
                                             int my_block_row, int my_block_col) const
{
  int lid_r[OPERATOR_MAX_NODES];
  int lid_c[OPERATOR_MAX_NODES];

  // ELEMENT: cell, DOFS: cell and node
  const std::vector<int>& node_row_inds = map.GhostIndices("node", my_block_row);
  const std::vector<int>& node_col_inds = map.GhostIndices("node", my_block_col);

  int ierr(0);
  AmanziMesh::Entity_ID_List nodes;
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    for (int n = 0; n != nnodes; ++n) {
      lid_r[n] = node_row_inds[nodes[n]];
      lid_c[n] = node_col_inds[nodes[n]];
    }
    ierr |= graph.InsertMyIndices(nnodes, lid_r, nnodes, lid_c);
  }
  ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble
* Apply the local matrices directly as schemas match.
****************************************************************** */
void Operator_Node::AssembleMatrixOp(const Op_Cell_Node& op,
                                     const SuperMap& map, MatrixFE& mat,
                                     int my_block_row, int my_block_col) const
{
  ASSERT(op.matrices.size() == ncells_owned);

  int lid_r[OPERATOR_MAX_NODES];
  int lid_c[OPERATOR_MAX_NODES];

  // ELEMENT: cell, DOFS: node and cell
  const std::vector<int>& node_row_inds = map.GhostIndices("node", my_block_row);
  const std::vector<int>& node_col_inds = map.GhostIndices("node", my_block_col);

  int ierr(0);
  AmanziMesh::Entity_ID_List nodes;
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    for (int n = 0; n != nnodes; ++n) {
      lid_r[n] = node_row_inds[nodes[n]];
      lid_c[n] = node_col_inds[nodes[n]];
    }

    ierr |= mat.SumIntoMyValues(lid_r, lid_c, op.matrices[c]);
  }
  ASSERT(!ierr);
}

}  // namespace Operators
}  // namespace Amanzi



