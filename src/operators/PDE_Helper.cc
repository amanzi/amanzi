/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Helper class for discrete PDE operators.
*/

#include "errors.hh"

#include "PDE_Helper.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Simple constructors.
****************************************************************** */
PDE_Helper::PDE_Helper(const Teuchos::RCP<Operator>& global_op) :
    global_op_(global_op)
{
  if (global_op == Teuchos::null) {
    Errors::Message msg("PDE_Helper: Constructor received null global operator");
    Exceptions::amanzi_throw(msg);
  }

  mesh_ = global_op_->Mesh();
  PopulateDimensions_();
}


PDE_Helper::PDE_Helper(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
     mesh_(mesh)
{
  PopulateDimensions_();
}


PDE_Helper::PDE_Helper(const Teuchos::RCP<AmanziMesh::Mesh>& mesh) :
     mesh_(mesh)
{
  PopulateDimensions_();
}


/* ******************************************************************
* Supporting private routines.
****************************************************************** */
void PDE_Helper::PopulateDimensions_()
{
  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);

  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::USED);

  if (mesh_->valid_edges()) {
    nedges_owned = mesh_->num_entities(AmanziMesh::EDGE, AmanziMesh::OWNED);
    nedges_wghost = mesh_->num_entities(AmanziMesh::EDGE, AmanziMesh::USED);
  }
}


/* ******************************************************************
* Apply BCs is typically called for the global Operator.
****************************************************************** */
void PDE_Helper::ApplyBCs_Face(const Teuchos::Ptr<BCs>& bcf, Teuchos::RCP<Op> op,
                               bool primary, bool eliminate)
{
  const std::vector<int>& bc_model = bcf->bc_model();
  const std::vector<double>& bc_value = bcf->bc_value();

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs, offset;

  CompositeVector& rhs = *global_op_->rhs();
  rhs.PutScalarGhosted(0.0);

  Teuchos::RCP<Epetra_MultiVector> rhs_face;
  if (primary) rhs_face = rhs.ViewComponent("face", true);

  const Schema& schema_row = global_op_->schema_row();
  const Schema& schema_col = global_op_->schema_col();

  for (int c = 0; c != ncells_owned; ++c) {
    WhetStone::DenseMatrix& Acell = op->matrices[c];
    int ncols = Acell.NumCols();
    int nrows = Acell.NumRows();

    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    // check for a boundary face
    bool found(false);
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      if (bc_model[f] == OPERATOR_BC_DIRICHLET) found = true;
    }
    if (!found) continue;

    // essential conditions for test functions
    schema_row.ComputeOffset(c, mesh_, offset);

    bool flag(true);
    int item(0);
    for (auto it = op->schema_row_.begin(); it != op->schema_row_.end(); ++it, ++item) {
      if (it->kind == AmanziMesh::FACE) {
        for (int n = 0; n != nfaces; ++n) {
          int f = faces[n];
          if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
            if (flag) {  // make a copy of elemental matrix
              op->matrices_shadow[c] = Acell;
              flag = false;
            }
            int noff(n + offset[item]);
            for (int m = 0; m < ncols; m++) Acell(noff, m) = 0.0;
          }
        }
      }
    }

    // essential zero conditions for trial functions
    schema_col.ComputeOffset(c, mesh_, offset);

    item = 0;
    for (auto it = op->schema_col_.begin(); it != op->schema_col_.end(); ++it, ++item) {
      if (it->kind == AmanziMesh::FACE) {
        for (int n = 0; n != nfaces; ++n) {
          int f = faces[n];
          double value = bc_value[f];

          if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
            if (flag) {  // make a copy of elemental matrix
              op->matrices_shadow[c] = Acell;
              flag = false;
            }

            int noff(n + offset[item]);
            WhetStone::DenseVector rhs_loc(nrows);

            if (eliminate) {
              for (int m = 0; m < nrows; m++) {
                rhs_loc(m) = -Acell(m, noff) * value;
                Acell(m, noff) = 0.0;
              }
            }

            if (primary) {
              rhs_loc(noff) = 0.0;
              (*rhs_face)[0][f] = value;
              Acell(noff, noff) = 1.0;
            }

            global_op_->AssembleVectorCellOp(c, schema_row, rhs_loc, rhs);
          }
        }
      }
    }
  } 

  rhs.GatherGhostedToMaster(Add);
}


/* ******************************************************************
* Apply BCs is typically called for the global Operator.
****************************************************************** */
void PDE_Helper::ApplyBCs_Node(const Teuchos::Ptr<BCs>& bcv, Teuchos::RCP<Op> op,
                               bool primary, bool eliminate)
{
  const std::vector<int>& bc_model = bcv->bc_model();
  const std::vector<AmanziGeometry::Point>& bc_value = bcv->bc_value_point();

  AmanziMesh::Entity_ID_List nodes, cells;
  std::vector<int> offset;

  CompositeVector& rhs = *global_op_->rhs();
  rhs.PutScalarGhosted(0.0);

  Teuchos::RCP<Epetra_MultiVector> rhs_node;
  if (primary) rhs_node = rhs.ViewComponent("node", true);

  int d = mesh_->space_dimension(); 
  const Schema& schema_row = global_op_->schema_row();
  const Schema& schema_col = global_op_->schema_col();

  for (int c = 0; c != ncells_owned; ++c) {
    WhetStone::DenseMatrix& Acell = op->matrices[c];
    int ncols = Acell.NumCols();
    int nrows = Acell.NumRows();

    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    // check for a boundary face
    bool found(false);
    for (int n = 0; n != nnodes; ++n) {
      int v = nodes[n];
      if (bc_model[v] == OPERATOR_BC_DIRICHLET) found = true;
    }
    if (!found) continue;

    // essential conditions for test functions
    op->schema_row_.ComputeOffset(c, mesh_, offset);

    bool flag(true);
    int item(0);
    for (auto it = op->schema_row_.begin(); it != op->schema_row_.end(); ++it, ++item) {
      if (it->kind == AmanziMesh::NODE) {
        for (int n = 0; n != nnodes; ++n) {
          int v = nodes[n];
          if (bc_model[v] == OPERATOR_BC_DIRICHLET) {
            if (flag) {  // make a copy of elemental matrix
              op->matrices_shadow[c] = Acell;
              flag = false;
            }
            for (int k = 0; k < d; ++k) {
              int noff(d*n + k + offset[item]);
              for (int m = 0; m < ncols; m++) Acell(noff, m) = 0.0;
            }
          }
        }
      }
    }

    // essential zero conditions for trial functions
    schema_col.ComputeOffset(c, mesh_, offset);

    item = 0;
    for (auto it = op->schema_col_.begin(); it != op->schema_col_.end(); ++it, ++item) {
      if (it->kind == AmanziMesh::NODE) {
        for (int n = 0; n != nnodes; ++n) {
          int v = nodes[n];
          AmanziGeometry::Point value = bc_value[v];

          if (bc_model[v] == OPERATOR_BC_DIRICHLET) {
            if (flag) {  // make a copy of elemental matrix
              op->matrices_shadow[c] = Acell;
              flag = false;
            }

            for (int k = 0; k < d; ++k) {
              int noff(d*n + k + offset[item]);
              WhetStone::DenseVector rhs_loc(nrows);

              if (eliminate) {
                for (int m = 0; m < nrows; m++) {
                  rhs_loc(m) = -Acell(m, noff) * value[k];
                  Acell(m, noff) = 0.0;
                }
              }

              if (primary) {
                mesh_->node_get_cells(v, AmanziMesh::USED, &cells);
                rhs_loc(noff) = 0.0;
                (*rhs_node)[k][v] = value[k];
                Acell(noff, noff) = 1.0 / cells.size();
              }

              global_op_->AssembleVectorCellOp(c, schema_row, rhs_loc, rhs);
            }
          }
        }
      }
    }
  } 

  rhs.GatherGhostedToMaster(Add);
}

}  // namespace Operators
}  // namespace Amanzi



