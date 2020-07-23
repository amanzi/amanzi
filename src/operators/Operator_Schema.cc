/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Operator those domain and range are defined by two schemas and 
  respected CVSs.
*/

#include "Epetra_Vector.h"

#include "DenseMatrix.hh"
#include "GraphFE.hh"
#include "MatrixFE.hh"
#include "MeshDefs.hh"
#include "SuperMap.hh"
#include "WhetStoneMeshUtils.hh"

#include "OperatorDefs.hh"
#include "OperatorUtils.hh"
#include "Operator_Schema.hh"
#include "Op_Cell_Schema.hh"
#include "Op_Face_Schema.hh"
#include "Op_Node_Schema.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Apply a source which may or may not have edge volume included already. 
****************************************************************** */
void Operator_Schema::UpdateRHS(const CompositeVector& source, bool volume_included)
{
  if (volume_included) {
    Operator::UpdateRHS(source);
  } else {
    AMANZI_ASSERT(false);
  }
}


/* ******************************************************************
* Visit methods for a matrix free mat-vec.
* Apply cell-based local matrices directly as schemas match.
****************************************************************** */
int Operator_Schema::ApplyMatrixFreeOp(const Op_Cell_Schema& op,
                                       const CompositeVector& X, CompositeVector& Y) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

  X.ScatterMasterToGhosted();
  Y.PutScalarGhosted(0.0);

  for (int c = 0; c != ncells_owned; ++c) {
    const WhetStone::DenseMatrix& A = op.matrices[c];
    int ncols = A.NumCols();
    int nrows = A.NumRows();
    WhetStone::DenseVector v(ncols), av(nrows);

    ExtractVectorCellOp(c, op.schema_col_, v, X);
    A.Multiply(v, av, false);
    AssembleVectorCellOp(c, op.schema_row_, av, Y);
  }

  Y.GatherGhostedToMaster(Add);
  return 0;
}


/* ******************************************************************
* Visit methods for a matrix free mat-vec.
* Apply face-based local matrices directly as schemas match.
****************************************************************** */
int Operator_Schema::ApplyMatrixFreeOp(const Op_Face_Schema& op,
                                       const CompositeVector& X, CompositeVector& Y) const
{
  AMANZI_ASSERT(op.matrices.size() == nfaces_owned);

  X.ScatterMasterToGhosted();
  Y.PutScalarGhosted(0.0);

  for (int f = 0; f != nfaces_owned; ++f) {
    const WhetStone::DenseMatrix& A = op.matrices[f];
    int ncols = A.NumCols();
    int nrows = A.NumRows();
    WhetStone::DenseVector v(ncols), av(nrows);

    ExtractVectorFaceOp(f, op.schema_col_, v, X);
    A.Multiply(v, av, false);
    AssembleVectorFaceOp(f, op.schema_row_, av, Y);
  }

  Y.GatherGhostedToMaster(Add);
  return 0;
}


/* ******************************************************************
* Visit methods for a matrix free mat-vec.
* Apply node-based local matrices directly as schemas match.
****************************************************************** */
int Operator_Schema::ApplyMatrixFreeOp(const Op_Node_Schema& op,
                                       const CompositeVector& X, CompositeVector& Y) const
{
  AMANZI_ASSERT(op.matrices.size() == nnodes_owned);

  X.ScatterMasterToGhosted();
  Y.PutScalarGhosted(0.0);

  for (int n = 0; n != nnodes_owned; ++n) {
    const WhetStone::DenseMatrix& A = op.matrices[n];
    int ncols = A.NumCols();
    int nrows = A.NumRows();
    WhetStone::DenseVector v(ncols), av(nrows);

    ExtractVectorNodeOp(n, op.schema_col_, v, X);
    A.Multiply(v, av, false);
    AssembleVectorNodeOp(n, op.schema_row_, av, Y);
  }

  Y.GatherGhostedToMaster(Add);
  return 0;
}


/* ******************************************************************
* Apply the local matrices directly as schemas match.
****************************************************************** */
int Operator_Schema::ApplyMatrixFreeOp(const Op_Node_Node& op,
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
* This method is mainly for debugging.
******************************************************************* */
int Operator_Schema::ApplyAssembled(const CompositeVector& X, CompositeVector& Y, double scalar) const
{
  X.ScatterMasterToGhosted();
  Y.PutScalarMasterAndGhosted(0.0);

  Epetra_Vector Xcopy(A_->RowMap());
  Epetra_Vector Ycopy(A_->RowMap());

  int ierr = CopyCompositeVectorToSuperVector(*smap_, X, Xcopy, schema_col_);
  ierr |= A_->Apply(Xcopy, Ycopy);
  ierr |= CopySuperVectorToCompositeVector(*smap_, Ycopy, Y, schema_row_);

  if (ierr) {
    Errors::Message msg;
    msg << "Operators: ApplyAssemble failed.\n";
    Exceptions::amanzi_throw(msg);
  }

  return ierr;
}



/* ******************************************************************
* Visit methods for symbolic assemble.
* Apply cell-based local matrices directly as schemas match.
****************************************************************** */
void Operator_Schema::SymbolicAssembleMatrixOp(const Op_Cell_Schema& op,
                                               const SuperMap& map, GraphFE& graph,
                                               int my_block_row, int my_block_col) const
{
  std::vector<int> lid_r, lid_c;
  AmanziMesh::Entity_ID_List entities;

  int num, ierr(0);
  AmanziMesh::Entity_kind kind;

  for (int c = 0; c != ncells_owned; ++c) {
    lid_c.clear();
    for (auto it = op.schema_col_.begin(); it != op.schema_col_.end(); ++it) {
      std::tie(kind, std::ignore, num) = *it;

      std::string name = schema_row_.KindToString(kind);
      WhetStone::cell_get_entities(*mesh_, c, kind, &entities);
      int nents = entities.size();
      AMANZI_ASSERT(nents > 0);

      for (int n = 0; n != nents; ++n) {
        int id = entities[n];
        for (int k = 0; k < num; ++k) {
          const std::vector<int>& col_inds = map.GhostIndices(my_block_col, name, k);
          lid_c.push_back(col_inds[id]);
        }
      }
    }

    lid_r.clear();
    for (auto it = op.schema_row_.begin(); it != op.schema_row_.end(); ++it) {
      std::tie(kind, std::ignore, num) = *it;

      std::string name = schema_row_.KindToString(kind);
      WhetStone::cell_get_entities(*mesh_, c, kind, &entities);
      int nents = entities.size();
      AMANZI_ASSERT(nents > 0);

      for (int n = 0; n != nents; ++n) {
        int id = entities[n];
        for (int k = 0; k < num; ++k) {
          const std::vector<int>& row_inds = map.GhostIndices(my_block_row, name, k);
          lid_r.push_back(row_inds[id]);
        }
      }
    }

    ierr |= graph.InsertMyIndices(lid_r.size(), lid_r.data(), lid_c.size(), lid_c.data());
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for symbolic assemble.
* Apply face-based local matrices directly as schemas match.
****************************************************************** */
void Operator_Schema::SymbolicAssembleMatrixOp(const Op_Face_Schema& op,
                                               const SuperMap& map, GraphFE& graph,
                                               int my_block_row, int my_block_col) const
{
  std::vector<int> lid_r, lid_c;
  AmanziMesh::Entity_ID_List cells;

  int ierr(0);
  for (int f = 0; f != nfaces_owned; ++f) {
    lid_r.clear();
    lid_c.clear();
    for (auto it = op.schema_col_.begin(); it != op.schema_col_.end(); ++it) {
      int num;
      AmanziMesh::Entity_kind kind;
      std::tie(kind, std::ignore, num) = *it;

      if (kind == AmanziMesh::CELL) {
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        int ncells = cells.size();

        for (int n = 0; n != ncells; ++n) {
          for (int k = 0; k < num; ++k) {
            const std::vector<int>& col_inds = map.GhostIndices(my_block_col, "cell", k);
            const std::vector<int>& row_inds = map.GhostIndices(my_block_row, "cell", k);

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
* Apply node-based local matrices directly as schemas match.
****************************************************************** */
void Operator_Schema::SymbolicAssembleMatrixOp(const Op_Node_Schema& op,
                                               const SuperMap& map, GraphFE& graph,
                                               int my_block_row, int my_block_col) const
{
  std::vector<int> lid_r, lid_c;
  AmanziMesh::Entity_ID_List cells;

  int ierr(0);
  for (int v = 0; v != nnodes_owned; ++v) {
    lid_r.clear();
    lid_c.clear();
    for (auto it = op.schema_col_.begin(); it != op.schema_col_.end(); ++it) {
      int num;
      AmanziMesh::Entity_kind kind;
      std::tie(kind, std::ignore, num) = *it;

      if (kind == AmanziMesh::CELL) {
        mesh_->node_get_cells(v, AmanziMesh::Parallel_type::ALL, &cells);
        int ncells = cells.size();

        for (int n = 0; n != ncells; ++n) {
          for (int k = 0; k < num; ++k) {
            const std::vector<int>& col_inds = map.GhostIndices(my_block_col, "cell", k);
            const std::vector<int>& row_inds = map.GhostIndices(my_block_row, "cell", k);

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
void Operator_Schema::SymbolicAssembleMatrixOp(const Op_Node_Node& op,
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
* Insert cell-based local matrices directly as schemas match.
****************************************************************** */
void Operator_Schema::AssembleMatrixOp(const Op_Cell_Schema& op,
                                       const SuperMap& map, MatrixFE& mat,
                                       int my_block_row, int my_block_col) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

  std::vector<int> lid_r, lid_c;
  AmanziMesh::Entity_ID_List entities;

  int num, ierr(0);
  AmanziMesh::Entity_kind kind;

  for (int c = 0; c != ncells_owned; ++c) {
    lid_c.clear();
    for (auto it = op.schema_col_.begin(); it != op.schema_col_.end(); ++it) {
      std::tie(kind, std::ignore, num) = *it;

      std::string name = schema_row_.KindToString(kind);
      WhetStone::cell_get_entities(*mesh_, c, kind, &entities);
      int nents = entities.size();
      AMANZI_ASSERT(nents > 0);

      for (int n = 0; n != nents; ++n) {
        int id = entities[n];
        for (int k = 0; k < num; ++k) {
          const std::vector<int>& col_inds = map.GhostIndices(my_block_col, name, k);
          lid_c.push_back(col_inds[id]);
        }
      }
    }

    lid_r.clear();
    for (auto it = op.schema_row_.begin(); it != op.schema_row_.end(); ++it) {
      std::tie(kind, std::ignore, num) = *it;

      std::string name = schema_row_.KindToString(kind);
      WhetStone::cell_get_entities(*mesh_, c, kind, &entities);
      int nents = entities.size();
      AMANZI_ASSERT(nents > 0);

      for (int n = 0; n != nents; ++n) {
        int id = entities[n];
        for (int k = 0; k < num; ++k) {
          const std::vector<int>& row_inds = map.GhostIndices(my_block_row, name, k);
          lid_r.push_back(row_inds[id]);
        }
      }
    }

    ierr |= mat.SumIntoMyValues(lid_r.data(), lid_c.data(), op.matrices[c]);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble
* Insert face-based local matrices directly as schemas match.
****************************************************************** */
void Operator_Schema::AssembleMatrixOp(const Op_Face_Schema& op,
                                       const SuperMap& map, MatrixFE& mat,
                                       int my_block_row, int my_block_col) const
{
  AMANZI_ASSERT(op.matrices.size() == nfaces_owned);

  std::vector<int> lid_r, lid_c;
  AmanziMesh::Entity_ID_List cells;

  int ierr(0);
  for (int f = 0; f != nfaces_owned; ++f) {
    lid_r.clear();
    lid_c.clear();
    for (auto it = op.schema_col_.begin(); it != op.schema_col_.end(); ++it) {
      int num;
      AmanziMesh::Entity_kind kind;
      std::tie(kind, std::ignore, num) = *it;

      if (kind == AmanziMesh::CELL) {
        mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
        int ncells = cells.size();

        for (int n = 0; n != ncells; ++n) {
          for (int k = 0; k < num; ++k) {
            const std::vector<int>& col_inds = map.GhostIndices(my_block_col, "cell", k);
            const std::vector<int>& row_inds = map.GhostIndices(my_block_row, "cell", k);

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
* Insert node-based local matrices directly as schemas match.
****************************************************************** */
void Operator_Schema::AssembleMatrixOp(const Op_Node_Schema& op,
                                       const SuperMap& map, MatrixFE& mat,
                                       int my_block_row, int my_block_col) const
{
  AMANZI_ASSERT(op.matrices.size() == nnodes_owned);

  std::vector<int> lid_r, lid_c;
  AmanziMesh::Entity_ID_List cells;

  int ierr(0);
  for (int v = 0; v != nnodes_owned; ++v) {
    lid_r.clear();
    lid_c.clear();
    for (auto it = op.schema_col_.begin(); it != op.schema_col_.end(); ++it) {
      int num;
      AmanziMesh::Entity_kind kind;
      std::tie(kind, std::ignore, num) = *it;

      if (kind == AmanziMesh::CELL) {
        mesh_->node_get_cells(v, AmanziMesh::Parallel_type::ALL, &cells);
        int ncells = cells.size();

        for (int n = 0; n != ncells; ++n) {
          for (int k = 0; k < num; ++k) {
            const std::vector<int>& col_inds = map.GhostIndices(my_block_col, "cell", k);
            const std::vector<int>& row_inds = map.GhostIndices(my_block_row, "cell", k);

            lid_c.push_back(col_inds[cells[n]]);
            lid_r.push_back(row_inds[cells[n]]);
          }
        }
      } else {
        AMANZI_ASSERT(false);
      }
    }

    ierr |= mat.SumIntoMyValues(lid_r.data(), lid_c.data(), op.matrices[v]);
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble
* Insert each diagonal values for nodes.
****************************************************************** */
void Operator_Schema::AssembleMatrixOp(const Op_Node_Node& op,
                                       const SuperMap& map, MatrixFE& mat,
                                       int my_block_row, int my_block_col) const
{
  const std::vector<int>& node_row_inds = map.GhostIndices(my_block_row, "node", 0);
  const std::vector<int>& node_col_inds = map.GhostIndices(my_block_col, "node", 0);

  int ierr(0);
  for (int v = 0; v != nnodes_owned; ++v) {
    int row = node_row_inds[v];
    int col = node_col_inds[v];

    for (int k = 0; k != op.diag->NumVectors(); ++k) {
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
void Operator_Schema::AssembleVectorCellOp(
    int c, const Schema& schema,
    const WhetStone::DenseVector& v, CompositeVector& X) const
{
  AmanziMesh::Entity_ID_List nodes, edges, faces;

  int m(0);
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    int num;
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, num) = *it;

    if (kind == AmanziMesh::NODE) {
      Epetra_MultiVector& Xn = *X.ViewComponent("node", true);

      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      for (int n = 0; n != nnodes; ++n) {
        for (int k = 0; k < num; ++k) {
          Xn[k][nodes[n]] += v(m++);
        }
      }
    }

    else if (kind == AmanziMesh::FACE) {
      Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      for (int n = 0; n != nfaces; ++n) {
        for (int k = 0; k < num; ++k) {
          Xf[k][faces[n]] += v(m++);
        }
      }
    }

    else if (kind == AmanziMesh::EDGE) {
      Epetra_MultiVector& Xe = *X.ViewComponent("edge", true);

      mesh_->cell_get_edges(c, &edges);
      int nedges = edges.size();

      for (int n = 0; n != nedges; ++n) {
        for (int k = 0; k < num; ++k) {
          Xe[k][edges[n]] += v(m++);
        }
      }
    }

    else if (kind == AmanziMesh::CELL) {
      Epetra_MultiVector& Xc = *X.ViewComponent("cell", true);

      for (int k = 0; k < num; ++k) {
        Xc[k][c] += v(m++);
      }
    }

    else {
      AMANZI_ASSERT(false);
    }
  }
}


/* ******************************************************************
* Assemble local vector to the global CV when the base is face.
****************************************************************** */
void Operator_Schema::AssembleVectorFaceOp(
    int f, const Schema& schema,
    const WhetStone::DenseVector& v, CompositeVector& X) const
{
  AmanziMesh::Entity_ID_List cells;

  int m(0);
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    int num;
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, num) = *it;

    if (kind == AmanziMesh::CELL) {
      Epetra_MultiVector& Xf = *X.ViewComponent("cell", true);

      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int ncells = cells.size();

      for (int n = 0; n != ncells; ++n) {
        for (int k = 0; k < num; ++k) {
          Xf[k][cells[n]] += v(m++);
        }
      }
    }
  }
}


/* ******************************************************************
* Assemble local vector to the global CV when the base is node.
****************************************************************** */
void Operator_Schema::AssembleVectorNodeOp(
    int n, const Schema& schema,
    const WhetStone::DenseVector& v, CompositeVector& X) const
{
  AmanziMesh::Entity_ID_List cells;

  int m(0);
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    int num;
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, num) = *it;

    if (kind == AmanziMesh::CELL) {
      Epetra_MultiVector& Xc = *X.ViewComponent("cell", true);

      mesh_->node_get_cells(n, AmanziMesh::Parallel_type::ALL, &cells);
      int ncells = cells.size();

      for (int i = 0; i != ncells; ++i) {
        for (int k = 0; k < num; ++k) {
          Xc[k][cells[i]] += v(m++);
        }
      }
    }
  }
}


/* ******************************************************************
* Extract local vector from the global CV when the base is cell.
****************************************************************** */
void Operator_Schema::ExtractVectorCellOp(
    int c, const Schema& schema,
    WhetStone::DenseVector& v, const CompositeVector& X) const
{
  AmanziMesh::Entity_ID_List nodes, edges, faces;

  int m(0);
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    int num;
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, num) = *it;

    if (kind == AmanziMesh::NODE) {
      const Epetra_MultiVector& Xn = *X.ViewComponent("node", true);

      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      for (int n = 0; n != nnodes; ++n) {
        for (int k = 0; k < num; ++k) {
          v(m++) = Xn[k][nodes[n]];
        }
      }
    }

    else if (kind == AmanziMesh::FACE) {
      const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

      mesh_->cell_get_faces(c, &faces);
      int nfaces = faces.size();

      for (int n = 0; n != nfaces; ++n) {
        for (int k = 0; k < num; ++k) {
          v(m++) = Xf[k][faces[n]];
        }
      }
    }

    else if (kind == AmanziMesh::EDGE) {
      const Epetra_MultiVector& Xe = *X.ViewComponent("edge", true);

      mesh_->cell_get_edges(c, &edges);
      int nedges = edges.size();

      for (int n = 0; n != nedges; ++n) {
        for (int k = 0; k < num; ++k) {
          v(m++) = Xe[k][edges[n]];
        }
      }
    }

    else if (kind == AmanziMesh::CELL) {
      const Epetra_MultiVector& Xc = *X.ViewComponent("cell", true);

      for (int k = 0; k < num; ++k) {
        v(m++) = Xc[k][c];
      }
    }

    else {
      AMANZI_ASSERT(false);
    }
  }
}


/* ******************************************************************
* Extract local vector from the global CV when the base is face.
****************************************************************** */
void Operator_Schema::ExtractVectorFaceOp(
    int f, const Schema& schema,
    WhetStone::DenseVector& v, const CompositeVector& X) const
{
  AmanziMesh::Entity_ID_List cells;

  int m(0);
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    int num;
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, num) = *it;

    if (kind == AmanziMesh::CELL) {
      const Epetra_MultiVector& Xf = *X.ViewComponent("cell", true);

      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      int ncells = cells.size();

      for (int n = 0; n != ncells; ++n) {
        for (int k = 0; k < num; ++k) {
          v(m++) = Xf[k][cells[n]];
        }
      }
    }
  }
}


/* ******************************************************************
* Extract local vector from the global CV when the base is node.
****************************************************************** */
void Operator_Schema::ExtractVectorNodeOp(
    int n, const Schema& schema,
    WhetStone::DenseVector& v, const CompositeVector& X) const
{
  AmanziMesh::Entity_ID_List cells;

  int m(0);
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    int num;
    AmanziMesh::Entity_kind kind;
    std::tie(kind, std::ignore, num) = *it;

    if (kind == AmanziMesh::CELL) {
      const Epetra_MultiVector& Xc = *X.ViewComponent("cell", true);

      mesh_->node_get_cells(n, AmanziMesh::Parallel_type::ALL, &cells);
      int ncells = cells.size();

      for (int i = 0; i != ncells; ++i) {
        for (int k = 0; k < num; ++k) {
          v(m++) = Xc[k][cells[i]];
        }
      }
    }
  }
}

}  // namespace Operators
}  // namespace Amanzi



