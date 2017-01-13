/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

// TPLs
#include "Epetra_Vector.h"

// Amanzi
#include "errors.hh"
#include "MatrixFE.hh"
#include "mfd3d_elasticity.hh"
#include "PreconditionerFactory.hh"
#include "WhetStoneDefs.hh"

// Amanzi::Operators
#include "Elasticity.hh"
#include "Op.hh"
#include "Op_Cell_Schema.hh"
#include "OperatorDefs.hh"
#include "Operator_Schema.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization of the operator, scalar coefficient.
****************************************************************** */
void Elasticity::SetTensorCoefficient(
    const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K)
{
  K_ = K;
}


/* ******************************************************************
* Calculate elemental matrices.
****************************************************************** */
void Elasticity::UpdateMatrices()
{
  WhetStone::MFD3D_Elasticity mfd(mesh_);
  AmanziMesh::Entity_ID_List nodes;
  int d = mesh_->space_dimension();

  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = 1.0;
  
  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();
    int nfaces = mesh_->cell_get_num_faces(c);
    int ndofs = d * nnodes + nfaces;

    WhetStone::DenseMatrix Acell(ndofs, ndofs);
    if (K_.get()) Kc = (*K_)[c];
    mfd.StiffnessMatrixNode2Face1(c, Kc, Acell);

    local_op_->matrices[c] = Acell;
  }
}


/* ******************************************************************
* Apply boundary conditions to the local matrices. We always zero-out
* matrix rows for essential test BCs. As to trial BCs, there are
* options: (a) eliminate or not, (b) if eliminate, then put 1 on
* the diagonal or not.
****************************************************************** */
void Elasticity::ApplyBCs(bool primary, bool eliminate)
{
  for (auto bc = bcs_trial_.begin(); bc != bcs_trial_.end(); ++bc) {
    if ((*bc)->type() == OPERATOR_BC_TYPE_FACE) {
      ApplyBCs_Face_(bc->ptr(), primary, eliminate);
    } else if ((*bc)->type() == OPERATOR_BC_TYPE_NODE) {
      ApplyBCs_Node_(bc->ptr(), primary, eliminate);
    }
  }
}


/* ******************************************************************
* Apply BCs on cell operators
****************************************************************** */
void Elasticity::ApplyBCs_Face_(const Teuchos::Ptr<BCs>& bcf,
                                bool primary, bool eliminate)
{
  const std::vector<int>& bc_model = bcf->bc_model();
  const std::vector<double>& bc_value = bcf->bc_value();

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs, offset;

  CompositeVector& rhs = *global_op_->rhs();
  Epetra_MultiVector& rhs_face = *rhs.ViewComponent("face", true);
  rhs.PutScalarGhosted(0.0);

  for (int c = 0; c != ncells_owned; ++c) {
    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];
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

    local_op_->schema_row_.ComputeOffset(c, mesh_, offset);

    // essential conditions for test functions
    int item(0);
    for (auto it = local_op_->schema_row_.begin(); it != local_op_->schema_row_.end(); ++it, ++item) {
      if (it->kind == AmanziMesh::FACE) {
        bool flag(true);
        for (int n = 0; n != nfaces; ++n) {
          int f = faces[n];
          if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
            if (flag) {  // make a copy of elemental matrix
              local_op_->matrices_shadow[c] = Acell;
              flag = false;
            }
            int noff(n + offset[item]);
            for (int m = 0; m < ncols; m++) Acell(noff, m) = 0.0;
          }
        }

        for (int n = 0; n != nfaces; ++n) {
          int f = faces[n];
          double value = bc_value[f];

          if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
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
              rhs_face[0][f] = value;
              Acell(noff, noff) = 1.0;
            }

            global_op_->AssembleVectorOp(c, global_schema_row_, rhs_loc, rhs);
          }
        }
      }
    }
  } 

  global_op_->rhs()->GatherGhostedToMaster("face", Add);
}


/* ******************************************************************
* Apply BCs on cell operators
****************************************************************** */
void Elasticity::ApplyBCs_Node_(const Teuchos::Ptr<BCs>& bcv,
                                bool primary, bool eliminate)
{
  const std::vector<int>& bc_model = bcv->bc_model();
  const std::vector<AmanziGeometry::Point>& bc_value = bcv->bc_value_point();

  AmanziMesh::Entity_ID_List nodes, cells;
  std::vector<int> offset;

  CompositeVector& rhs = *global_op_->rhs();
  Epetra_MultiVector& rhs_node = *rhs.ViewComponent("node", true);
  rhs.PutScalarGhosted(0.0);

  int d = mesh_->space_dimension(); 

  for (int c = 0; c != ncells_owned; ++c) {
    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];
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

    local_op_->schema_row_.ComputeOffset(c, mesh_, offset);

    // essential conditions for test functions
    int item(0);
    for (auto it = local_op_->schema_row_.begin(); it != local_op_->schema_row_.end(); ++it, ++item) {
      if (it->kind == AmanziMesh::NODE) {
        bool flag(true);
        for (int n = 0; n != nnodes; ++n) {
          int v = nodes[n];
          if (bc_model[v] == OPERATOR_BC_DIRICHLET) {
            if (flag) {  // make a copy of elemental matrix
              local_op_->matrices_shadow[c] = Acell;
              flag = false;
            }
            for (int k = 0; k < d; ++k) {
              int noff(d*n + k + offset[item]);
              for (int m = 0; m < ncols; m++) Acell(noff, m) = 0.0;
            }
          }
        }

        for (int n = 0; n != nnodes; ++n) {
          int v = nodes[n];
          AmanziGeometry::Point value = bc_value[v];

          if (bc_model[v] == OPERATOR_BC_DIRICHLET) {
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
                rhs_node[k][v] = value[k];
                Acell(noff, noff) = 1.0 / cells.size();
              }

              global_op_->AssembleVectorOp(c, global_schema_row_, rhs_loc, rhs);
            }
          }
        }
      }
    }
  } 

  global_op_->rhs()->GatherGhostedToMaster("node", Add);
}


/* ******************************************************************
* Put here stuff that has to be done in constructor.
****************************************************************** */
void Elasticity::InitElasticity_(Teuchos::ParameterList& plist)
{
  // Read schema for the mimetic discretization method.
  Teuchos::ParameterList& schema_list = plist.sublist("schema");

  std::vector<std::string> name;
  if (schema_list.isParameter("location")) {
    name = schema_list.get<Teuchos::Array<std::string> >("location").toVector();
  } else {
    Errors::Message msg;
    msg << "Elasticity: schema->location is missing.";
    Exceptions::amanzi_throw(msg);
  }

  std::vector<std::string> type;
  if (schema_list.isParameter("type")) {
    type = schema_list.get<Teuchos::Array<std::string> >("type").toVector();
  } else {
    Errors::Message msg;
    msg << "Elasticity: schema->type is missing.";
    Exceptions::amanzi_throw(msg);
  }

  std::vector<int> ndofs;
  if (schema_list.isParameter("number")) {
    ndofs = schema_list.get<Teuchos::Array<int> >("number").toVector();
  } else {
    Errors::Message msg;
    msg << "Elasticity: schema->number is missing.";
    Exceptions::amanzi_throw(msg);
  }

  // Create schema and save it.
  Schema my_schema;
  my_schema.SetBase(OPERATOR_SCHEMA_BASE_CELL);
  for (int i = 0; i < name.size(); i++) {
    my_schema.AddItem(my_schema.StringToKind(name[i]), SCHEMA_DOFS_SCALAR, ndofs[i]);
  }
  my_schema.Finalize(mesh_);

  // create or check the existing Operator
  local_schema_row_ = my_schema;
  local_schema_row_ = my_schema;

  if (global_op_ == Teuchos::null) {
    global_schema_col_ = my_schema;
    global_schema_row_ = my_schema;

    // build the CVS from the global schema
    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true);

    for (int i = 0; i < name.size(); i++) {
      cvs->AddComponent(name[i], my_schema.StringToKind(name[i]), ndofs[i]);
    }

    global_op_ = Teuchos::rcp(new Operator_Schema(cvs, cvs, plist, my_schema, my_schema));
  } else {
    // constructor was given an Operator
    global_schema_col_ = global_op_->schema_col();
    global_schema_row_ = global_op_->schema_row();
    mesh_ = global_op_->DomainMap().Mesh();
  }

  // create the local Op and register it with the global Operator
  local_op_ = Teuchos::rcp(new Op_Cell_Schema(my_schema, my_schema, mesh_));
  global_op_->OpPushBack(local_op_);
  
  // mesh info
  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);

  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  K_ = Teuchos::null;
}

}  // namespace Operators
}  // namespace Amanzi
