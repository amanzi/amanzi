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
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

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
  const Epetra_MultiVector& Xn = *X.ViewComponent("node");
  Epetra_MultiVector& Yn = *Y.ViewComponent("node");

  for (int i = 0; i < Xn.NumVectors(); ++i) {
    for (int v = 0; v != nnodes_owned; ++v) {
      Yn[i][v] += Xn[i][v] * (*op.diag)[i][v];
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
  std::vector<int> lid_r(cell_max_nodes);
  std::vector<int> lid_c(cell_max_nodes);

  // ELEMENT: cell, DOFS: cell and node
  const std::vector<int>& node_row_inds = map.GhostIndices(my_block_row, "node", 0);
  const std::vector<int>& node_col_inds = map.GhostIndices(my_block_col, "node", 0);

  int ierr(0);
  AmanziMesh::Entity_ID_List nodes;
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    for (int n = 0; n != nnodes; ++n) {
      lid_r[n] = node_row_inds[nodes[n]];
      lid_c[n] = node_col_inds[nodes[n]];
    }
    ierr |= graph.InsertMyIndices(nnodes, lid_r.data(), nnodes, lid_c.data());
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for symbolic assemble.
* Insert the diagonal at nodes
****************************************************************** */
void Operator_Node::SymbolicAssembleMatrixOp(const Op_Node_Node& op,
                                             const SuperMap& map, GraphFE& graph,
                                             int my_block_row, int my_block_col) const
{
  const std::vector<int>& node_row_inds = map.GhostIndices(my_block_row, "node", 0);
  const std::vector<int>& node_col_inds = map.GhostIndices(my_block_col, "node", 0);

  int ierr(0);
  for (int v = 0; v != nnodes_owned; ++v) {
    int row = node_row_inds[v];
    int col = node_col_inds[v];

    ierr |= graph.InsertMyIndices(row, 1, &col);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble
* Apply the local matrices directly as schemas match.
****************************************************************** */
void Operator_Node::AssembleMatrixOp(const Op_Cell_Node& op,
                                     const SuperMap& map, MatrixFE& mat,
                                     int my_block_row, int my_block_col) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

  std::vector<int> lid_r(cell_max_nodes);
  std::vector<int> lid_c(cell_max_nodes);

  // ELEMENT: cell, DOFS: node and cell
  const std::vector<int>& node_row_inds = map.GhostIndices(my_block_row, "node", 0);
  const std::vector<int>& node_col_inds = map.GhostIndices(my_block_col, "node", 0);

  int ierr(0);
  AmanziMesh::Entity_ID_List nodes;
  for (int c = 0; c != ncells_owned; ++c) {
    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    for (int n = 0; n != nnodes; ++n) {
      lid_r[n] = node_row_inds[nodes[n]];
      lid_c[n] = node_col_inds[nodes[n]];
    }

    ierr |= mat.SumIntoMyValues(lid_r.data(), lid_c.data(), op.matrices[c]);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble
* Insert each diagonal values for edges.
****************************************************************** */
void Operator_Node::AssembleMatrixOp(const Op_Node_Node& op,
                                     const SuperMap& map, MatrixFE& mat,
                                     int my_block_row, int my_block_col) const
{
  AMANZI_ASSERT(op.diag->NumVectors() == 1);

  const std::vector<int>& node_row_inds = map.GhostIndices(my_block_row, "node", 0);
  const std::vector<int>& node_col_inds = map.GhostIndices(my_block_col, "node", 0);

  int ierr(0);
  for (int v = 0; v != nnodes_owned; ++v) {
    int row = node_row_inds[v];
    int col = node_col_inds[v];

    ierr |= mat.SumIntoMyValues(row, 1, &(*op.diag)[0][v], &col);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Operator> Operator_Node::Clone() const {
  return Teuchos::rcp(new Operator_Node(*this));
}

}  // namespace Operators
}  // namespace Amanzi



