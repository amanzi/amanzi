/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "Epetra_Vector.h"

#include "DenseMatrix.hh"
#include "GraphFE.hh"
#include "MatrixFE.hh"
#include "MeshDefs.hh"
#include "PreconditionerFactory.hh"
#include "SuperMap.hh"

#include "OperatorDefs.hh"
#include "OperatorUtils.hh"
#include "Operator_Schema.hh"
#include "Op_Cell_Schema.hh"
#include "Op_Face_Schema.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
 * Apply a source which may or may not have edge volume included already.
 ****************************************************************** */
void
Operator_Schema::UpdateRHS(const CompositeVector& source, bool volume_included)
{
  if (volume_included) {
    Operator::UpdateRHS(source);
  } else {
    AMANZI_ASSERT(false);
  }
}


/* ******************************************************************
 * Apply the local matrices directly as schemas match.
 ****************************************************************** */
int
Operator_Schema::ApplyMatrixFreeOp(const Op_Cell_Schema& op,
                                   const CompositeVector& X,
                                   CompositeVector& Y) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

  X.ScatterMasterToGhosted();
  Y.putScalarGhosted(0.0);

  for (int c = 0; c != ncells_owned; ++c) {
    const WhetStone::DenseMatrix& Acell = op.matrices[c];
    int ncols = Acell.NumCols();
    int nrows = Acell.NumRows();
    WhetStone::DenseVector v(ncols), av(nrows);

    ExtractVectorCellOp(c, op.schema_col_, v, X);
    Acell.elementWiseMultiply(v, av, false);
    AssembleVectorCellOp(c, op.schema_row_, av, Y);
  }

  Y.GatherGhostedToMaster(Add);
  return 0;
}


int
Operator_Schema::ApplyMatrixFreeOp(const Op_Face_Schema& op,
                                   const CompositeVector& X,
                                   CompositeVector& Y) const
{
  AMANZI_ASSERT(op.matrices.size() == nfaces_owned);

  X.ScatterMasterToGhosted();
  Y.putScalarGhosted(0.0);

  for (int f = 0; f != nfaces_owned; ++f) {
    const WhetStone::DenseMatrix& Aface = op.matrices[f];
    int ncols = Aface.NumCols();
    int nrows = Aface.NumRows();
    WhetStone::DenseVector v(ncols), av(nrows);

    ExtractVectorFaceOp(f, op.schema_col_, v, X);
    Aface.elementWiseMultiply(v, av, false);
    AssembleVectorFaceOp(f, op.schema_row_, av, Y);
  }

  Y.GatherGhostedToMaster(Add);
  return 0;
}


/* ******************************************************************
 * Apply the local matrices directly as schemas match.
 ****************************************************************** */
int
Operator_Schema::ApplyMatrixFreeOp(const Op_Node_Node& op,
                                   const CompositeVector& X,
                                   CompositeVector& Y) const
{
  const Epetra_MultiVector& Xn = *X.ViewComponent("node");
  Epetra_MultiVector& Yn = *Y.ViewComponent("node");

  for (int i = 0; i < Xn.getNumVectors(); ++i) {
    for (int v = 0; v != nnodes_owned; ++v) {
      Yn[i][v] += Xn[i][v] * (*op.diag)[i][v];
    }
  }
  return 0;
}


/* ******************************************************************
 * Apply the local matrices directly as schemas match.
 ****************************************************************** */
int
Operator_Schema::ApplyTransposeMatrixFreeOp(const Op_Cell_Schema& op,
                                            const CompositeVector& X,
                                            CompositeVector& Y) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

  X.ScatterMasterToGhosted();
  Y.putScalarGhosted(0.0);

  for (int c = 0; c != ncells_owned; ++c) {
    const WhetStone::DenseMatrix& Acell = op.matrices[c];
    int ncols = Acell.NumCols();
    int nrows = Acell.NumRows();
    WhetStone::DenseVector v(nrows), av(ncols);

    ExtractVectorCellOp(c, op.schema_row_, v, X);
    Acell.elementWiseMultiply(v, av, true);
    AssembleVectorCellOp(c, op.schema_col_, av, Y);
  }

  Y.GatherGhostedToMaster(Add);
  return 0;
}


/* ******************************************************************
 * Parallel preconditioner: Y = P * X.
 ******************************************************************* */
int
Operator_Schema::ApplyInverse(const CompositeVector& X,
                              CompositeVector& Y) const
{
  Epetra_Vector Xcopy(A_->RowMap());
  Epetra_Vector Ycopy(A_->RowMap());

  int ierr = CopyCompositeVectorToSuperVector(*smap_, X, Xcopy, schema_col_);
  ierr |= preconditioner_->applyInverse(Xcopy, Ycopy);
  ierr |= CopySuperVectorToCompositeVector(*smap_, Ycopy, Y, schema_row_);

  return ierr;
}


/* ******************************************************************
 * This method is mainly for debugging.
 ******************************************************************* */
int
Operator_Schema::ApplyAssembled(const CompositeVector& X, CompositeVector& Y,
                                double scalar) const
{
  X.ScatterMasterToGhosted();
  Y.putScalarMasterAndGhosted(0.0);

  Epetra_Vector Xcopy(A_->RowMap());
  Epetra_Vector Ycopy(A_->RowMap());

  int ierr = CopyCompositeVectorToSuperVector(*smap_, X, Xcopy, schema_col_);
  ierr |= A_->apply(Xcopy, Ycopy);
  ierr |= CopySuperVectorToCompositeVector(*smap_, Ycopy, Y, schema_row_);

  if (ierr) {
    Errors::Message msg;
    msg << "Operators: ApplyAssemble failed.\n";
    Exceptions::amanzi_throw(msg);
  }

  return ierr;
}


/* ******************************************************************
 * Create structure of a global square matrix.
 ****************************************************************** */
void
Operator_Schema::SymbolicAssembleMatrix()
{
  // Create the supermap given a space (set of possible schemas) and a
  // specific schema (assumed/checked to be consistent with the space).
  smap_ = CreateSuperMap(*cvs_col_, schema_col_);

  // create the graph
  int row_size = MaxRowSize(*mesh_, schema_col_);
  Teuchos::RCP<GraphFE> graph = Teuchos::rcp(new GraphFE(
    smap_->getMap(), smap_->GhostedMap(), smap_->GhostedMap(), row_size));

  // fill the graph
  Operator::SymbolicAssembleMatrix(*smap_, *graph, 0, 0);

  // Completing and optimizing the graphs
  int ierr = graph->FillComplete(smap_->getMap(), smap_->getMap());
  AMANZI_ASSERT(!ierr);

  // create global matrix
  Amat_ = Teuchos::rcp(new MatrixFE(graph));
  A_ = Amat_->Matrix();
}


/* ******************************************************************
 * Visit methods for symbolic assemble.
 * Apply the local matrices directly as schemas match.
 ****************************************************************** */
void
Operator_Schema::SymbolicAssembleMatrixOp(const Op_Cell_Schema& op,
                                          const SuperMap& map, GraphFE& graph,
                                          int my_block_row,
                                          int my_block_col) const
{
  std::vector<int> lid_r, lid_c;
  AmanziMesh::Entity_ID_List nodes, edges, faces;

  int ierr(0);
  for (int c = 0; c != ncells_owned; ++c) {
    lid_r.clear();
    lid_c.clear();
    for (auto it = op.schema_col_.begin(); it != op.schema_col_.end(); ++it) {
      if (it->kind == AmanziMesh::NODE) {
        mesh_->cell_get_nodes(c, &nodes);
        int nnodes = nodes.size();

        for (int n = 0; n != nnodes; ++n) {
          int v = nodes[n];
          for (int k = 0; k < it->num; ++k) {
            const std::vector<int>& col_inds =
              map.GhostIndices(my_block_col, "node", k);
            const std::vector<int>& row_inds =
              map.GhostIndices(my_block_row, "node", k);

            lid_c.push_back(col_inds[v]);
            lid_r.push_back(row_inds[v]);
          }
        }
      }

      else if (it->kind == AmanziMesh::FACE) {
        mesh_->cell_get_faces(c, &faces);
        int nfaces = faces.size();

        for (int n = 0; n != nfaces; ++n) {
          int f = faces[n];
          for (int k = 0; k < it->num; ++k) {
            const std::vector<int>& col_inds =
              map.GhostIndices(my_block_col, "face", k);
            const std::vector<int>& row_inds =
              map.GhostIndices(my_block_row, "face", k);

            lid_c.push_back(col_inds[f]);
            lid_r.push_back(row_inds[f]);
          }
        }
      }

      else if (it->kind == AmanziMesh::EDGE) {
        mesh_->cell_get_edges(c, &edges);
        int nedges = edges.size();

        for (int n = 0; n != nedges; ++n) {
          int e = edges[n];
          for (int k = 0; k < it->num; ++k) {
            const std::vector<int>& col_inds =
              map.GhostIndices(my_block_col, "edge", k);
            const std::vector<int>& row_inds =
              map.GhostIndices(my_block_row, "edge", k);

            lid_c.push_back(col_inds[e]);
            lid_r.push_back(row_inds[e]);
          }
        }
      }

      else if (it->kind == AmanziMesh::CELL) {
        for (int k = 0; k < it->num; ++k) {
          const std::vector<int>& col_inds =
            map.GhostIndices(my_block_col, "cell", k);
          const std::vector<int>& row_inds =
            map.GhostIndices(my_block_row, "cell", k);

          lid_c.push_back(col_inds[c]);
          lid_r.push_back(row_inds[c]);
        }
      }

      else {
        AMANZI_ASSERT(false);
      }
    }

    int m = lid_c.size();
    ierr |= graph.InsertMyIndices(m, lid_r.data(), m, lid_c.data());
  }
  AMANZI_ASSERT(!ierr);
}


void
Operator_Schema::SymbolicAssembleMatrixOp(const Op_Face_Schema& op,
                                          const SuperMap& map, GraphFE& graph,
                                          int my_block_row,
                                          int my_block_col) const
{
  std::vector<int> lid_r, lid_c;
  AmanziMesh::Entity_ID_List cells;

  int ierr(0);
  for (int f = 0; f != nfaces_owned; ++f) {
    lid_r.clear();
    lid_c.clear();
    for (auto it = op.schema_col_.begin(); it != op.schema_col_.end(); ++it) {
      if (it->kind == AmanziMesh::CELL) {
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        int ncells = cells.size();

        for (int n = 0; n != ncells; ++n) {
          for (int k = 0; k < it->num; ++k) {
            const std::vector<int>& col_inds =
              map.GhostIndices(my_block_col, "cell", k);
            const std::vector<int>& row_inds =
              map.GhostIndices(my_block_row, "cell", k);

            lid_c.push_back(col_inds[cells[n]]);
            lid_r.push_back(row_inds[cells[n]]);
          }
        }
      } else {
        AMANZI_ASSERT(false);
      }
    }

    int m = lid_c.size();
    ierr |= graph.InsertMyIndices(m, lid_r.data(), m, lid_c.data());
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
 * Visit methods for symbolic assemble.
 * Insert a diagonal matrix on nodes.
 ****************************************************************** */
void
Operator_Schema::SymbolicAssembleMatrixOp(const Op_Node_Node& op,
                                          const SuperMap& map, GraphFE& graph,
                                          int my_block_row,
                                          int my_block_col) const
{
  const std::vector<int>& node_row_inds =
    map.GhostIndices(my_block_row, "node", 0);
  const std::vector<int>& node_col_inds =
    map.GhostIndices(my_block_col, "node", 0);

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
void
Operator_Schema::AssembleMatrixOp(const Op_Cell_Schema& op, const SuperMap& map,
                                  MatrixFE& mat, int my_block_row,
                                  int my_block_col) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

  std::vector<int> lid_r, lid_c;
  AmanziMesh::Entity_ID_List nodes, edges, faces;

  int ierr(0);
  for (int c = 0; c != ncells_owned; ++c) {
    lid_r.clear();
    lid_c.clear();
    for (auto it = op.schema_col_.begin(); it != op.schema_col_.end(); ++it) {
      if (it->kind == AmanziMesh::NODE) {
        mesh_->cell_get_nodes(c, &nodes);
        int nnodes = nodes.size();

        for (int n = 0; n != nnodes; ++n) {
          int v = nodes[n];
          for (int k = 0; k < it->num; ++k) {
            const std::vector<int>& col_inds =
              map.GhostIndices(my_block_col, "node", k);
            const std::vector<int>& row_inds =
              map.GhostIndices(my_block_row, "node", k);

            lid_c.push_back(col_inds[nodes[n]]);
            lid_r.push_back(row_inds[nodes[n]]);
          }
        }
      }

      else if (it->kind == AmanziMesh::FACE) {
        mesh_->cell_get_faces(c, &faces);
        int nfaces = faces.size();

        for (int n = 0; n != nfaces; ++n) {
          int f = faces[n];
          for (int k = 0; k < it->num; ++k) {
            const std::vector<int>& col_inds =
              map.GhostIndices(my_block_col, "face", k);
            const std::vector<int>& row_inds =
              map.GhostIndices(my_block_row, "face", k);

            lid_c.push_back(col_inds[f]);
            lid_r.push_back(row_inds[f]);
          }
        }
      }

      else if (it->kind == AmanziMesh::EDGE) {
        mesh_->cell_get_edges(c, &edges);
        int nedges = edges.size();

        for (int n = 0; n != nedges; ++n) {
          int e = edges[n];
          for (int k = 0; k < it->num; ++k) {
            const std::vector<int>& col_inds =
              map.GhostIndices(my_block_col, "edge", k);
            const std::vector<int>& row_inds =
              map.GhostIndices(my_block_row, "edge", k);

            lid_c.push_back(col_inds[e]);
            lid_r.push_back(row_inds[e]);
          }
        }
      }

      else if (it->kind == AmanziMesh::CELL) {
        for (int k = 0; k < it->num; ++k) {
          const std::vector<int>& col_inds =
            map.GhostIndices(my_block_col, "cell", k);
          const std::vector<int>& row_inds =
            map.GhostIndices(my_block_row, "cell", k);

          lid_c.push_back(col_inds[c]);
          lid_r.push_back(row_inds[c]);
        }
      }

      else {
        AMANZI_ASSERT(false);
      }
    }

    ierr |= mat.SumIntoMyValues(lid_r.data(), lid_c.data(), op.matrices[c]);
  }
  AMANZI_ASSERT(!ierr);
}


void
Operator_Schema::AssembleMatrixOp(const Op_Face_Schema& op, const SuperMap& map,
                                  MatrixFE& mat, int my_block_row,
                                  int my_block_col) const
{
  AMANZI_ASSERT(op.matrices.size() == nfaces_owned);

  std::vector<int> lid_r, lid_c;
  AmanziMesh::Entity_ID_List cells;

  int ierr(0);
  for (int f = 0; f != nfaces_owned; ++f) {
    lid_r.clear();
    lid_c.clear();
    for (auto it = op.schema_col_.begin(); it != op.schema_col_.end(); ++it) {
      if (it->kind == AmanziMesh::CELL) {
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        int ncells = cells.size();

        for (int n = 0; n != ncells; ++n) {
          for (int k = 0; k < it->num; ++k) {
            const std::vector<int>& col_inds =
              map.GhostIndices(my_block_col, "cell", k);
            const std::vector<int>& row_inds =
              map.GhostIndices(my_block_row, "cell", k);

            lid_c.push_back(col_inds[cells[n]]);
            lid_r.push_back(row_inds[cells[n]]);
          }
        }
      } else {
        AMANZI_ASSERT(false);
      }
    }

    ierr |= mat.SumIntoMyValues(lid_r.data(), lid_c.data(), op.matrices[f]);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
 * Visit methods for assemble
 * Insert each diagonal values for nodes.
 ****************************************************************** */
void
Operator_Schema::AssembleMatrixOp(const Op_Node_Node& op, const SuperMap& map,
                                  MatrixFE& mat, int my_block_row,
                                  int my_block_col) const
{
  const std::vector<int>& node_row_inds =
    map.GhostIndices(my_block_row, "node", 0);
  const std::vector<int>& node_col_inds =
    map.GhostIndices(my_block_col, "node", 0);

  int ierr(0);
  for (int v = 0; v != nnodes_owned; ++v) {
    int row = node_row_inds[v];
    int col = node_col_inds[v];

    for (int k = 0; k != op.diag->getNumVectors(); ++k) {
      ierr |= mat.SumIntoMyValues(row, 1, &(*op.diag)[k][v], &col);
      row++;
      col++;
    }
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
 * Assemble local vector to the global CV when the base is cell.
 ****************************************************************** */
void
Operator_Schema::AssembleVectorCellOp(int c, const Schema& schema,
                                      const WhetStone::DenseVector& v,
                                      CompositeVector& X) const
{
  AmanziMesh::Entity_ID_List nodes, edges, faces;

  int m(0);
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    if (it->kind == AmanziMesh::NODE) {
      Epetra_MultiVector& Xn = *X.ViewComponent("node", true);

      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      for (int n = 0; n != nnodes; ++n) {
        for (int k = 0; k < it->num; ++k) { Xn[k][nodes[n]] += v(m++); }
      }
    }

    else if (it->kind == AmanziMesh::FACE) {
      Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      for (int n = 0; n != nfaces; ++n) {
        for (int k = 0; k < it->num; ++k) { Xf[k][faces[n]] += v(m++); }
      }
    }

    else if (it->kind == AmanziMesh::EDGE) {
      Epetra_MultiVector& Xe = *X.ViewComponent("edge", true);

      mesh_->cell_get_edges(c, &edges);
      int nedges = edges.size();

      for (int n = 0; n != nedges; ++n) {
        for (int k = 0; k < it->num; ++k) { Xe[k][edges[n]] += v(m++); }
      }
    }

    else if (it->kind == AmanziMesh::CELL) {
      Epetra_MultiVector& Xc = *X.ViewComponent("cell", true);

      for (int k = 0; k < it->num; ++k) { Xc[k][c] += v(m++); }
    }

    else {
      AMANZI_ASSERT(false);
    }
  }
}


/* ******************************************************************
 * Assemble local vector to the global CV when the base is face.
 ****************************************************************** */
void
Operator_Schema::AssembleVectorFaceOp(int f, const Schema& schema,
                                      const WhetStone::DenseVector& v,
                                      CompositeVector& X) const
{
  AmanziMesh::Entity_ID_List cells;

  int m(0);
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    if (it->kind == AmanziMesh::CELL) {
      Epetra_MultiVector& Xf = *X.ViewComponent("cell", true);

      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int ncells = cells.size();

      for (int n = 0; n != ncells; ++n) {
        for (int k = 0; k < it->num; ++k) { Xf[k][cells[n]] += v(m++); }
      }
    }
  }
}


/* ******************************************************************
 * Extract local vector from the global CV when the base is cell.
 ****************************************************************** */
void
Operator_Schema::ExtractVectorCellOp(int c, const Schema& schema,
                                     WhetStone::DenseVector& v,
                                     const CompositeVector& X) const
{
  AmanziMesh::Entity_ID_List nodes, edges, faces;

  int m(0);
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    if (it->kind == AmanziMesh::NODE) {
      const Epetra_MultiVector& Xn = *X.ViewComponent("node", true);

      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      for (int n = 0; n != nnodes; ++n) {
        for (int k = 0; k < it->num; ++k) { v(m++) = Xn[k][nodes[n]]; }
      }
    }

    else if (it->kind == AmanziMesh::FACE) {
      const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      for (int n = 0; n != nfaces; ++n) {
        for (int k = 0; k < it->num; ++k) { v(m++) = Xf[k][faces[n]]; }
      }
    }

    else if (it->kind == AmanziMesh::EDGE) {
      const Epetra_MultiVector& Xe = *X.ViewComponent("edge", true);

      mesh_->cell_get_edges(c, &edges);
      int nedges = edges.size();

      for (int n = 0; n != nedges; ++n) {
        for (int k = 0; k < it->num; ++k) { v(m++) = Xe[k][edges[n]]; }
      }
    }

    else if (it->kind == AmanziMesh::CELL) {
      const Epetra_MultiVector& Xc = *X.ViewComponent("cell", true);

      for (int k = 0; k < it->num; ++k) { v(m++) = Xc[k][c]; }
    }

    else {
      AMANZI_ASSERT(false);
    }
  }
}


/* ******************************************************************
 * Extract local vector from the global CV when the base is face.
 ****************************************************************** */
void
Operator_Schema::ExtractVectorFaceOp(int f, const Schema& schema,
                                     WhetStone::DenseVector& v,
                                     const CompositeVector& X) const
{
  AmanziMesh::Entity_ID_List cells;

  int m(0);
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    if (it->kind == AmanziMesh::CELL) {
      const Epetra_MultiVector& Xf = *X.ViewComponent("cell", true);

      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int ncells = cells.size();

      for (int n = 0; n != ncells; ++n) {
        for (int k = 0; k < it->num; ++k) { v(m++) = Xf[k][cells[n]]; }
      }
    }
  }
}

} // namespace Operators
} // namespace Amanzi
