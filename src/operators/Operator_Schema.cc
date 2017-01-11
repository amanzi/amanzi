/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  Operator those range is defined by two schemas and respected CVSs.
*/

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
    ASSERT(false);
  }
}


/* ******************************************************************
* Apply the local matrices directly as schemas match.
****************************************************************** */
int Operator_Schema::ApplyMatrixFreeOp(const Op_Cell_Schema& op,
                                       const CompositeVector& X, CompositeVector& Y) const
{
  ASSERT(op.matrices.size() == ncells_owned);

  X.ScatterMasterToGhosted();
  Y.PutScalarGhosted(0.0);

  for (int c = 0; c != ncells_owned; ++c) {
    const WhetStone::DenseMatrix& Acell = op.matrices[c];
    int ncols = Acell.NumCols();
    int nrows = Acell.NumRows();
    WhetStone::DenseVector v(ncols), av(nrows);

    ExtractVectorOp(c, op.schema_col_, v, X);
    Acell.Multiply(v, av, false);
    AssembleVectorOp(c, op.schema_row_, av, Y);
  }

  Y.GatherGhostedToMaster(Add);
  return 0;
}


/* ******************************************************************
* Parallel preconditioner: Y = P * X.
******************************************************************* */
int Operator_Schema::ApplyInverse(const CompositeVector& X, CompositeVector& Y) const
{
  Y = X;
  return 0;
  Epetra_Vector Xcopy(A_->RowMap());
  Epetra_Vector Ycopy(A_->RowMap());
  int ierr = CopyCompositeVectorToSuperVector(*smap_, X, Xcopy, schema_col_);

  ierr |= preconditioner_->ApplyInverse(Xcopy, Ycopy);
  ierr |= CopySuperVectorToCompositeVector(*smap_, Ycopy, Y, schema_row_);
  ASSERT(!ierr);

  return ierr;
}


/* ******************************************************************
* Create structure of a global square matrix.
****************************************************************** */
void Operator_Schema::SymbolicAssembleMatrix()
{
  // Create the supermap given a space (set of possible schemas) and a
  // specific schema (assumed/checked to be consistent with the space).
  smap_ = CreateSuperMap(*cvs_col_, schema_col_);

  // create the graph
  int row_size = MaxRowSize(*mesh_, schema_col_);
  Teuchos::RCP<GraphFE> graph = Teuchos::rcp(new GraphFE(smap_->Map(),
          smap_->GhostedMap(), smap_->GhostedMap(), row_size));

  // fill the graph
  Operator::SymbolicAssembleMatrix(*smap_, *graph, 0, 0);

  // Completing and optimizing the graphs
  int ierr = graph->FillComplete(smap_->Map(), smap_->Map());
  ASSERT(!ierr);

  // create global matrix
  Amat_ = Teuchos::rcp(new MatrixFE(graph));
  A_ = Amat_->Matrix();
}


/* ******************************************************************
* Visit methods for symbolic assemble.
* Apply the local matrices directly as schemas match.
****************************************************************** */
void Operator_Schema::SymbolicAssembleMatrixOp(const Op_Cell_Schema& op,
                                               const SuperMap& map, GraphFE& graph,
                                               int my_block_row, int my_block_col) const
{
  int lid_c[OPERATOR_MAX_EDGES];
  int lid_r[OPERATOR_MAX_EDGES];

  AmanziMesh::Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  int ierr(0);
  for (int c = 0; c != ncells_owned; ++c) {
    int m(0);
    for (auto it = op.schema_col_.begin(); it != op.schema_col_.end(); ++it) {
      if (it->kind == AmanziMesh::NODE) {
        mesh_->cell_get_nodes(c, &nodes);
        int nnodes = nodes.size();

        for (int k = 0; k < it->num; ++k) {
          const std::vector<int>& col_inds = map.GhostIndices("node", k);
          const std::vector<int>& row_inds = map.GhostIndices("node", k);

          for (int n = 0; n != nnodes; ++n) {
            lid_c[m] = col_inds[nodes[n]];
            lid_r[m] = row_inds[nodes[n]];
            m++;
          }
        }
      }

      if (it->kind == AmanziMesh::FACE) {
        mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
        int nfaces = faces.size();

        for (int k = 0; k < it->num; ++k) {
          const std::vector<int>& col_inds = map.GhostIndices("face", k);
          const std::vector<int>& row_inds = map.GhostIndices("face", k);

          for (int n = 0; n != nfaces; ++n) {
            lid_c[m] = col_inds[faces[n]];
            lid_r[m] = row_inds[faces[n]];
            m++;
          }
        }
      }
    }

    ierr |= graph.InsertMyIndices(m, lid_r, m, lid_c);
  }
  ASSERT(!ierr);
}


/* ******************************************************************
* Visit methods for assemble
* Apply the local matrices directly as schemas match.
****************************************************************** */
void Operator_Schema::AssembleMatrixOp(const Op_Cell_Schema& op,
                                       const SuperMap& map, MatrixFE& mat,
                                       int my_block_row, int my_block_col) const
{
  ASSERT(op.matrices.size() == ncells_owned);

  int lid_r[OPERATOR_MAX_EDGES];
  int lid_c[OPERATOR_MAX_EDGES];

  AmanziMesh::Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  int ierr(0);
  for (int c = 0; c != ncells_owned; ++c) {
    int m(0);
    for (auto it = op.schema_col_.begin(); it != op.schema_col_.end(); ++it) {
      if (it->kind == AmanziMesh::NODE) {
        mesh_->cell_get_nodes(c, &nodes);
        int nnodes = nodes.size();

        for (int k = 0; k < it->num; ++k) {
          const std::vector<int>& col_inds = map.GhostIndices("node", k);
          const std::vector<int>& row_inds = map.GhostIndices("node", k);

          for (int n = 0; n != nnodes; ++n) {
            lid_c[m] = col_inds[nodes[n]];
            lid_r[m] = row_inds[nodes[n]];
            m++;
          }
        }
      }

      if (it->kind == AmanziMesh::FACE) {
        mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
        int nfaces = faces.size();

        for (int k = 0; k < it->num; ++k) {
          const std::vector<int>& col_inds = map.GhostIndices("face", k);
          const std::vector<int>& row_inds = map.GhostIndices("face", k);

          for (int n = 0; n != nfaces; ++n) {
            lid_c[m] = col_inds[faces[n]];
            lid_r[m] = row_inds[faces[n]];
            m++;
          }
        }
      }
    }

    ierr |= mat.SumIntoMyValues(lid_r, lid_c, op.matrices[c]);
  }
  ASSERT(!ierr);
}


/* ******************************************************************
* Assemble local vector to the global CV.
****************************************************************** */
void Operator_Schema::AssembleVectorOp(
    int c, const Schema& schema,
    const WhetStone::DenseVector& v, CompositeVector& X) const
{
  AmanziMesh::Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  int m(0);
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    if (it->kind == AmanziMesh::NODE) {
      Epetra_MultiVector& Xn = *X.ViewComponent("node", true);

      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      for (int n = 0; n != nnodes; ++n) {
        for (int k = 0; k < it->num; ++k) {
          Xn[k][nodes[n]] += v(m++);
        }
      }
    }

    if (it->kind == AmanziMesh::FACE) {
      Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();

      for (int n = 0; n != nfaces; ++n) {
        for (int k = 0; k < it->num; ++k) {
          Xf[k][faces[n]] += v(m++);
        }
      }
    }
  }
}


/* ******************************************************************
* Extract local vector to the global CV.
****************************************************************** */
void Operator_Schema::ExtractVectorOp(
    int c, const Schema& schema,
    WhetStone::DenseVector& v, const CompositeVector& X) const
{
  AmanziMesh::Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  int m(0);
  for (auto it = schema.begin(); it != schema.end(); ++it) {
    if (it->kind == AmanziMesh::NODE) {
      const Epetra_MultiVector& Xn = *X.ViewComponent("node", true);

      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      for (int n = 0; n != nnodes; ++n) {
        for (int k = 0; k < it->num; ++k) {
          v(m++) = Xn[k][nodes[n]];
        }
      }
    }

    if (it->kind == AmanziMesh::FACE) {
      const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

      mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
      int nfaces = faces.size();

      for (int n = 0; n != nfaces; ++n) {
        for (int k = 0; k < it->num; ++k) {
          v(m++) = Xf[k][faces[n]];
        }
      }
    }
  }
}

}  // namespace Operators
}  // namespace Amanzi



