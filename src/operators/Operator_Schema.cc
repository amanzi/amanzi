/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Operator those domain and range are defined by two schemas and
  respected CVSs.
*/

#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"

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
#include "Op_MeshInjection.hh"
#include "Op_Node_Node.hh"
#include "SchemaUtils.hh"

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
* Visit methods for a matrix free mat-vec.
* Apply cell-based local matrices directly as schemas match.
****************************************************************** */
int
Operator_Schema::ApplyMatrixFreeOp(const Op_Cell_Schema& op,
                                   const CompositeVector& X,
                                   CompositeVector& Y) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);
  for (int c = 0; c != ncells_owned; ++c) {
    const WhetStone::DenseMatrix& A = op.matrices[c];
    int ncols = A.NumCols();
    int nrows = A.NumRows();
    WhetStone::DenseVector v(ncols), av(nrows);

    ExtractVectorCellOp(c, *mesh_, op.schema_col(), v, X);
    A.Multiply(v, av, false);
    AssembleVectorCellOp(c, *mesh_, op.schema_row(), av, Y);
  }
  return 0;
}


/* ******************************************************************
* Visit methods for a matrix free mat-vec.
* Apply face-based local matrices directly as schemas match.
****************************************************************** */
int
Operator_Schema::ApplyMatrixFreeOp(const Op_Face_Schema& op,
                                   const CompositeVector& X,
                                   CompositeVector& Y) const
{
  AMANZI_ASSERT(op.matrices.size() == nfaces_owned);

  for (int f = 0; f != nfaces_owned; ++f) {
    const WhetStone::DenseMatrix& A = op.matrices[f];
    int ncols = A.NumCols();
    int nrows = A.NumRows();
    WhetStone::DenseVector v(ncols), av(nrows);

    ExtractVectorFaceOp(f, *mesh_, op.schema_col(), v, X);
    A.Multiply(v, av, false);
    AssembleVectorFaceOp(f, *mesh_, op.schema_row(), av, Y);
  }
  return 0;
}


/* ******************************************************************
* Visit methods for a matrix free mat-vec.
* Apply node-based local matrices directly as schemas match.
****************************************************************** */
int
Operator_Schema::ApplyMatrixFreeOp(const Op_Node_Schema& op,
                                   const CompositeVector& X,
                                   CompositeVector& Y) const
{
  AMANZI_ASSERT(op.matrices.size() == nnodes_owned);

  for (int n = 0; n != nnodes_owned; ++n) {
    const WhetStone::DenseMatrix& A = op.matrices[n];
    int ncols = A.NumCols();
    int nrows = A.NumRows();
    WhetStone::DenseVector v(ncols), av(nrows);

    ExtractVectorNodeOp(n, *mesh_, op.schema_col(), v, X);
    A.Multiply(v, av, false);
    AssembleVectorNodeOp(n, *mesh_, op.schema_row(), av, Y);
  }

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

  for (int i = 0; i < Xn.NumVectors(); ++i) {
    for (int v = 0; v != nnodes_owned; ++v) { Yn[i][v] += Xn[i][v] * (*op.diag)[i][v]; }
  }
  return 0;
}


/* ******************************************************************
* Apply the local matrices directly.
****************************************************************** */
int
Operator_Schema::ApplyMatrixFreeOp(const Op_MeshInjection& op,
                                   const CompositeVector& X,
                                   CompositeVector& Y) const
{
  auto col_comp_name = AmanziMesh::entity_kind_string(std::get<0>(*op.schema_col().begin()));
  const Epetra_MultiVector& X_vec = *X.ViewComponent(col_comp_name, false);

  auto row_comp_name = AmanziMesh::entity_kind_string(std::get<0>(*op.schema_row().begin()));
  Epetra_MultiVector& Y_vec = *Y.ViewComponent(row_comp_name, false);

  if (!op.transpose) {
    // restrict/inject from X into Y -- Y += A*X, where A is short and long
    Epetra_MultiVector X_restricted(Y_vec);
    Epetra_Import restriction(*op.injection, X_vec.Map());
    int ierr = X_restricted.Import(X_vec, restriction, Insert);
    AMANZI_ASSERT(!ierr);
    Y_vec.Multiply(1.0, *op.diag, X_restricted, 1.0);
  } else {
    // prolongate from X into Y -- Y += A*X, where A is tall and skinny
    Epetra_MultiVector Y_restricted(X_vec);
    Y_restricted.PutScalar(0.);
    Y_restricted.Multiply(1.0, *op.diag, X_vec, 0.0);

    Epetra_Import prolongation(Y_vec.Map(), *op.injection);
    int ierr = Y_vec.Import(Y_restricted, prolongation, Epetra_AddLocalAlso);
    AMANZI_ASSERT(!ierr);
  }
  return 0;
}


/* ******************************************************************
* This method is mainly for debugging.
******************************************************************* */
int
Operator_Schema::ApplyAssembled(const CompositeVector& X, CompositeVector& Y, double scalar) const
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
void
Operator_Schema::SymbolicAssembleMatrixOp(const Op_Cell_Schema& op,
                                          const SuperMap& map,
                                          GraphFE& graph,
                                          int my_block_row,
                                          int my_block_col) const
{
  std::vector<int> lid_r, lid_c;
  AmanziMesh::Entity_ID_List entities;

  int num, ierr(0);
  AmanziMesh::Entity_kind kind;

  for (int c = 0; c != ncells_owned; ++c) {
    lid_c.clear();
    for (auto it = op.schema_col().begin(); it != op.schema_col().end(); ++it) {
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
    for (auto it = op.schema_row().begin(); it != op.schema_row().end(); ++it) {
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
void
Operator_Schema::SymbolicAssembleMatrixOp(const Op_Face_Schema& op,
                                          const SuperMap& map,
                                          GraphFE& graph,
                                          int my_block_row,
                                          int my_block_col) const
{
  std::vector<int> lid_r, lid_c;
  AmanziMesh::Entity_ID_List cells;

  int ierr(0);
  for (int f = 0; f != nfaces_owned; ++f) {
    lid_r.clear();
    lid_c.clear();
    for (auto it = op.schema_col().begin(); it != op.schema_col().end(); ++it) {
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
void
Operator_Schema::SymbolicAssembleMatrixOp(const Op_Node_Schema& op,
                                          const SuperMap& map,
                                          GraphFE& graph,
                                          int my_block_row,
                                          int my_block_col) const
{
  std::vector<int> lid_r, lid_c;
  AmanziMesh::Entity_ID_List cells;

  int ierr(0);
  for (int v = 0; v != nnodes_owned; ++v) {
    lid_r.clear();
    lid_c.clear();
    for (auto it = op.schema_col().begin(); it != op.schema_col().end(); ++it) {
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
void
Operator_Schema::SymbolicAssembleMatrixOp(const Op_Node_Node& op,
                                          const SuperMap& map,
                                          GraphFE& graph,
                                          int my_block_row,
                                          int my_block_col) const
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
* Symbolic assemble an injection
****************************************************************** */
void
Operator_Schema::SymbolicAssembleMatrixOp(const Op_MeshInjection& op,
                                          const SuperMap& map,
                                          GraphFE& graph,
                                          int my_block_row,
                                          int my_block_col) const
{
  auto row_entity_kind = std::get<0>(*op.schema_row().begin());
  auto col_entity_kind = std::get<0>(*op.schema_col().begin());
  auto row_comp_name = AmanziMesh::entity_kind_string(row_entity_kind);
  auto col_comp_name = AmanziMesh::entity_kind_string(col_entity_kind);
  const std::vector<int>& row_inds = map.GhostIndices(my_block_row, row_comp_name, 0);
  const std::vector<int>& col_inds = map.GhostIndices(my_block_col, col_comp_name, 0);

  auto row_entity_map = op.get_row_mesh().map(row_entity_kind, false);
  auto col_entity_map = op.get_col_mesh().map(col_entity_kind, false);

  int ierr(0);
  if (!op.transpose) {
    AMANZI_ASSERT(row_entity_map.NumMyElements() == op.injection->NumMyElements());
    for (int row_lid = 0; row_lid != row_entity_map.NumMyElements(); ++row_lid) {
      auto col_lid = col_entity_map.LID(op.injection->GID(row_lid));
      int col = col_inds[col_lid];
      ierr |= graph.InsertMyIndices(row_inds[row_lid], 1, &col);
    }
  } else {
    AMANZI_ASSERT(col_entity_map.NumMyElements() == op.injection->NumMyElements());
    for (int col_lid = 0; col_lid != col_entity_map.NumMyElements(); ++col_lid) {
      auto row_lid = row_entity_map.LID(op.injection->GID(col_lid));
      int col = col_inds[col_lid];
      ierr |= graph.InsertMyIndices(row_inds[row_lid], 1, &col);
    }
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble
* Insert cell-based local matrices directly as schemas match.
****************************************************************** */
void
Operator_Schema::AssembleMatrixOp(const Op_Cell_Schema& op,
                                  const SuperMap& map,
                                  MatrixFE& mat,
                                  int my_block_row,
                                  int my_block_col) const
{
  AMANZI_ASSERT(op.matrices.size() == ncells_owned);

  std::vector<int> lid_r, lid_c;
  AmanziMesh::Entity_ID_List entities;

  int num, ierr(0);
  AmanziMesh::Entity_kind kind;

  for (int c = 0; c != ncells_owned; ++c) {
    lid_c.clear();
    for (auto it = op.schema_col().begin(); it != op.schema_col().end(); ++it) {
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
    for (auto it = op.schema_row().begin(); it != op.schema_row().end(); ++it) {
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
void
Operator_Schema::AssembleMatrixOp(const Op_Face_Schema& op,
                                  const SuperMap& map,
                                  MatrixFE& mat,
                                  int my_block_row,
                                  int my_block_col) const
{
  AMANZI_ASSERT(op.matrices.size() == nfaces_owned);

  std::vector<int> lid_r, lid_c;
  AmanziMesh::Entity_ID_List cells;

  int ierr(0);
  for (int f = 0; f != nfaces_owned; ++f) {
    lid_r.clear();
    lid_c.clear();
    for (auto it = op.schema_col().begin(); it != op.schema_col().end(); ++it) {
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
void
Operator_Schema::AssembleMatrixOp(const Op_Node_Schema& op,
                                  const SuperMap& map,
                                  MatrixFE& mat,
                                  int my_block_row,
                                  int my_block_col) const
{
  AMANZI_ASSERT(op.matrices.size() == nnodes_owned);

  std::vector<int> lid_r, lid_c;
  AmanziMesh::Entity_ID_List cells;

  int ierr(0);
  for (int v = 0; v != nnodes_owned; ++v) {
    lid_r.clear();
    lid_c.clear();
    for (auto it = op.schema_col().begin(); it != op.schema_col().end(); ++it) {
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
void
Operator_Schema::AssembleMatrixOp(const Op_Node_Node& op,
                                  const SuperMap& map,
                                  MatrixFE& mat,
                                  int my_block_row,
                                  int my_block_col) const
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
* Assemble an injection
****************************************************************** */
void
Operator_Schema::AssembleMatrixOp(const Op_MeshInjection& op,
                                  const SuperMap& map,
                                  MatrixFE& mat,
                                  int my_block_row,
                                  int my_block_col) const
{
  auto row_entity_kind = std::get<0>(*op.schema_row().begin());
  auto col_entity_kind = std::get<0>(*op.schema_col().begin());
  auto row_comp_name = AmanziMesh::entity_kind_string(row_entity_kind);
  auto col_comp_name = AmanziMesh::entity_kind_string(col_entity_kind);
  const std::vector<int>& row_inds = map.GhostIndices(my_block_row, row_comp_name, 0);
  const std::vector<int>& col_inds = map.GhostIndices(my_block_col, col_comp_name, 0);

  auto row_entity_map = op.get_row_mesh().map(row_entity_kind, false);
  auto col_entity_map = op.get_col_mesh().map(col_entity_kind, false);

  int ierr(0);
  if (!op.transpose) {
    AMANZI_ASSERT(row_entity_map.NumMyElements() == op.injection->NumMyElements());
    for (int row_lid = 0; row_lid != row_entity_map.NumMyElements(); ++row_lid) {
      auto col_lid = col_entity_map.LID(op.injection->GID(row_lid));
      ierr |=
        mat.SumIntoMyValues(row_inds[row_lid], 1, &(*op.diag)[0][row_lid], &col_inds[col_lid]);
    }
  } else {
    AMANZI_ASSERT(col_entity_map.NumMyElements() == op.injection->NumMyElements());
    for (int col_lid = 0; col_lid != col_entity_map.NumMyElements(); ++col_lid) {
      auto row_lid = row_entity_map.LID(op.injection->GID(col_lid));
      ierr |=
        mat.SumIntoMyValues(row_inds[row_lid], 1, &(*op.diag)[0][col_lid], &col_inds[col_lid]);
    }
  }
  AMANZI_ASSERT(!ierr);
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Operator>
Operator_Schema::Clone() const
{
  return Teuchos::rcp(new Operator_Schema(*this));
}

} // namespace Operators
} // namespace Amanzi
