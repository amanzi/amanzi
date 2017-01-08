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

#include "DenseMatrix.hh"
#include "Op_Edge_Edge.hh"
#include "Op_Cell_Edge.hh"

#include "SuperMap.hh"
#include "GraphFE.hh"
#include "MatrixFE.hh"

#include "OperatorDefs.hh"
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

  AmanziMesh::Entity_ID_List nodes, faces;
  std::vector<int> dirs;

  for (int c = 0; c != ncells_owned; ++c) {
    const WhetStone::DenseMatrix& Acell = op.matrices[c];
    int ncols = Acell.NumCols();
    int nrows = Acell.NumRows();
    WhetStone::DenseVector v(ncols), av(nrows);

    // extract local vector
    int m(0);
    for (auto it = op.schema_col_.begin(); it != op.schema_col_.end(); ++it) {
      if (it->location == OPERATOR_SCHEMA_DOFS_NODE) {
        const Epetra_MultiVector& Xn = *X.ViewComponent("node", true);

        mesh_->cell_get_nodes(c, &nodes);
        int nnodes = nodes.size();

        for (int n = 0; n != nnodes; ++n) {
          for (int k = 0; k < it->num; ++k) {
            v(m++) = Xn[k][nodes[n]];
          }
        }
      }

      if (it->location == OPERATOR_SCHEMA_DOFS_FACE) {
        const Epetra_MultiVector& Xf = *X.ViewComponent("face", true);

        mesh_->cell_get_nodes(c, &faces);
        int nfaces = faces.size();

        for (int n = 0; n != nfaces; ++n) {
          for (int k = 0; k < it->num; ++k) {
            v(m++) = Xf[k][faces[n]];
          }
        }
      }
    }

    Acell.Multiply(v, av, false);

    // assemble the local vector into the global one
    m = 0;
    for (auto it = op.schema_row_.begin(); it != op.schema_row_.end(); ++it) {
      if (it->location == OPERATOR_SCHEMA_DOFS_NODE) {
        Epetra_MultiVector& Yn = *Y.ViewComponent("node", true);

        mesh_->cell_get_nodes(c, &nodes);
        int nnodes = nodes.size();

        for (int n = 0; n != nnodes; ++n) {
          for (int k = 0; k < it->num; ++k) {
            Yn[k][nodes[n]] += av(m++);
          }
        }
      }

      if (it->location == OPERATOR_SCHEMA_DOFS_FACE) {
        Epetra_MultiVector& Yf = *Y.ViewComponent("face", true);

        mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
        int nfaces = faces.size();

        for (int n = 0; n != nfaces; ++n) {
          for (int k = 0; k < it->num; ++k) {
            Yf[k][faces[n]] += av(m++);
          }
        }
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
      if (it->location == OPERATOR_SCHEMA_DOFS_NODE) {
        mesh_->cell_get_nodes(c, &nodes);
        int nnodes = nodes.size();

        for (int k = 0; k < it->num; ++k) {
          const std::vector<int>& col_inds = map.GhostIndices("node", k);
          const std::vector<int>& row_inds = map.GhostIndices("node", k);

          for (int n = 0; n != nnodes; ++n) {
            lid_c[m++] = col_inds[nodes[n]];
            lid_r[m++] = row_inds[nodes[n]];
          }
        }
      }

      if (it->location == OPERATOR_SCHEMA_DOFS_FACE) {
        mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
        int nfaces = faces.size();

        for (int k = 0; k < it->num; ++k) {
          const std::vector<int>& col_inds = map.GhostIndices("face", k);
          const std::vector<int>& row_inds = map.GhostIndices("face", k);

          for (int n = 0; n != nfaces; ++n) {
            lid_c[m++] = col_inds[faces[n]];
            lid_r[m++] = row_inds[faces[n]];
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
      if (it->location == OPERATOR_SCHEMA_DOFS_NODE) {
        mesh_->cell_get_nodes(c, &nodes);
        int nnodes = nodes.size();

        for (int k = 0; k < it->num; ++k) {
          const std::vector<int>& col_inds = map.GhostIndices("node", k);
          const std::vector<int>& row_inds = map.GhostIndices("node", k);

          for (int n = 0; n != nnodes; ++n) {
            lid_c[m++] = col_inds[nodes[n]];
            lid_r[m++] = row_inds[nodes[n]];
          }
        }
      }

      if (it->location == OPERATOR_SCHEMA_DOFS_FACE) {
        mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
        int nfaces = faces.size();

        for (int k = 0; k < it->num; ++k) {
          const std::vector<int>& col_inds = map.GhostIndices("face", k);
          const std::vector<int>& row_inds = map.GhostIndices("face", k);

          for (int n = 0; n != nfaces; ++n) {
            lid_c[m++] = col_inds[faces[n]];
            lid_r[m++] = row_inds[faces[n]];
          }
        }
      }
    }

    ierr |= mat.SumIntoMyValues(lid_r, lid_c, op.matrices[c]);
  }
  ASSERT(!ierr);
}

}  // namespace Operators
}  // namespace Amanzi



