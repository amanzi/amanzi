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

#include "WhetStoneDefs.hh"

#include "ParallelCommunication.hh"
#include "PDE_HelperDiscretization.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Simple constructors.
****************************************************************** */
PDE_HelperDiscretization::PDE_HelperDiscretization(const Teuchos::RCP<Operator>& global_op) :
    global_op_(global_op)
{
  if (global_op == Teuchos::null) {
    Errors::Message msg("PDE_HelperDiscretization: Constructor received null global operator");
    Exceptions::amanzi_throw(msg);
  }

  mesh_ = global_op_->Mesh();
  PopulateDimensions_();
}


PDE_HelperDiscretization::PDE_HelperDiscretization(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
     mesh_(mesh)
{
  PopulateDimensions_();
}


PDE_HelperDiscretization::PDE_HelperDiscretization(const Teuchos::RCP<AmanziMesh::Mesh>& mesh) :
     mesh_(mesh)
{
  PopulateDimensions_();
}


/* ******************************************************************
* Supporting private routines.
****************************************************************** */
void PDE_HelperDiscretization::PopulateDimensions_()
{
  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);

  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);

  if (mesh_->valid_edges()) {
    nedges_owned = mesh_->num_entities(AmanziMesh::EDGE, AmanziMesh::Parallel_type::OWNED);
    nedges_wghost = mesh_->num_entities(AmanziMesh::EDGE, AmanziMesh::Parallel_type::ALL);
  }
}


/* ******************************************************************
* Replace container of local matrices with another container.
****************************************************************** */
void PDE_HelperDiscretization::set_local_op(const Teuchos::RCP<Op>& op)
{
  if (global_op_.get()) {
    if (local_op_.get()) {
      auto index = std::find(global_op_->begin(), global_op_->end(), local_op_) - global_op_->begin();
      if (index != global_op_->size()) {
        global_op_->OpPushBack(op);
      } else {
        global_op_->OpReplace(op, index);
      }
    } else {
      global_op_->OpPushBack(op);
    }
  }
  local_op_ = op;
}


/* ******************************************************************
* Apply boundary conditions to the local matrices. 
* NOTE: We always zero-out matrix rows for essential test BCs. 
****************************************************************** */
void PDE_HelperDiscretization::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  auto base = global_op_->schema_row().base();

  for (auto bc : bcs_trial_) {
    bool missing(true);

    if (bc->type() == WhetStone::DOF_Type::SCALAR ||
        bc->type() == WhetStone::DOF_Type::NORMAL_COMPONENT ||
        bc->type() == WhetStone::DOF_Type::MOMENT) {
      if (base == AmanziMesh::CELL) {
        ApplyBCs_Cell_Scalar_(*bc, local_op_, primary, eliminate, essential_eqn);
        missing = false;
      }
    } 
    else if (bc->type() == WhetStone::DOF_Type::POINT) {
      if (base == AmanziMesh::CELL) {
        ApplyBCs_Cell_Point_(*bc, local_op_, primary, eliminate, essential_eqn);
        missing = false;
      } else if (base == AmanziMesh::NODE) {
        missing = false;  // for testing
      }
    }
    else if (bc->type() == WhetStone::DOF_Type::VECTOR) {
      if (base == AmanziMesh::CELL) {
        ApplyBCs_Cell_Vector_(*bc, local_op_, primary, eliminate, essential_eqn);
        missing = false;
      }
    }

    if (missing) {
      Errors::Message msg("PDE_HelperDiscretization: Unsupported boundary condition.\n");
      Exceptions::amanzi_throw(msg);
    }
  }
}


/* ******************************************************************
* Apply BCs of scalar type.
****************************************************************** */
void PDE_HelperDiscretization::ApplyBCs_Cell_Scalar_(
    const BCs& bc, Teuchos::RCP<Op> op,
    bool primary, bool eliminate, bool essential_eqn)
{
  const std::vector<int>& bc_model = bc.bc_model();
  const std::vector<double>& bc_value = bc.bc_value();

  AmanziMesh::Entity_ID_List entities, cells;
  std::vector<int> offset;

  CompositeVector& rhs = *global_op_->rhs();
  rhs.PutScalarGhosted(0.0);

  const Schema& schema_row = global_op_->schema_row();
  const Schema& schema_col = global_op_->schema_col();

  AmanziMesh::Entity_kind kind = bc.kind();
  Teuchos::RCP<Epetra_MultiVector> rhs_kind;
  if (primary) rhs_kind = rhs.ViewComponent(schema_row.KindToString(kind), true);

  for (int c = 0; c != ncells_owned; ++c) {
    WhetStone::DenseMatrix& Acell = op->matrices[c];
    int ncols = Acell.NumCols();
    int nrows = Acell.NumRows();

    int nents_owned(0);
    if (kind == AmanziMesh::FACE) {
      mesh_->cell_get_faces(c, &entities);
      nents_owned = ncells_owned;
    } else if (kind == AmanziMesh::EDGE) {
      mesh_->cell_get_edges(c, &entities);
      nents_owned = nedges_owned;
    } else if (kind == AmanziMesh::NODE) {
      mesh_->cell_get_nodes(c, &entities);
      nents_owned = nnodes_owned;
    }
    int nents = entities.size();

    // check for a boundary entity
    bool found(false);
    for (int n = 0; n != nents; ++n) {
      int x = entities[n];
      if (bc_model[x] == OPERATOR_BC_DIRICHLET) found = true;
    }
    if (!found) continue;

    // essential conditions for test functions
    schema_row.ComputeOffset(c, mesh_, offset);

    bool flag(true);
    int item(0);
    AmanziMesh::Entity_kind itkind;

    for (auto it = op->schema_row_.begin(); it != op->schema_row_.end(); ++it, ++item) {
      std::tie(itkind, std::ignore, std::ignore) = *it;

      if (itkind == kind) {
        for (int n = 0; n != nents; ++n) {
          int x = entities[n];
          if (bc_model[x] == OPERATOR_BC_DIRICHLET) {
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
      std::tie(itkind, std::ignore, std::ignore) = *it;

      if (itkind == kind) {
        for (int n = 0; n != nents; ++n) {
          int x = entities[n];
          double value = bc_value[x];

          if (bc_model[x] == OPERATOR_BC_DIRICHLET) {
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

            if (essential_eqn) {
              rhs_loc(noff) = 0.0;
              (*rhs_kind)[0][x] = (x < nents_owned) ? value : 0.0;

              if (kind == AmanziMesh::FACE) {
                mesh_->face_get_cells(x, AmanziMesh::Parallel_type::ALL, &cells);
              } else if (kind == AmanziMesh::NODE) {
                mesh_->node_get_cells(x, AmanziMesh::Parallel_type::ALL, &cells);
              } else if (kind == AmanziMesh::EDGE) {
                mesh_->edge_get_cells(x, AmanziMesh::Parallel_type::ALL, &cells);
              }
              Acell(noff, noff) = 1.0 / cells.size();
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
* Apply BCs of point type. The code is limited to node DOFs.
****************************************************************** */
void PDE_HelperDiscretization::ApplyBCs_Cell_Point_(
    const BCs& bc, Teuchos::RCP<Op> op,
    bool primary, bool eliminate, bool essential_eqn)
{
  const std::vector<int>& bc_model = bc.bc_model();
  const std::vector<AmanziGeometry::Point>& bc_value = bc.bc_value_point();

  AmanziMesh::Entity_ID_List nodes, cells;
  std::vector<int> offset;

  CompositeVector& rhs = *global_op_->rhs();
  rhs.PutScalarGhosted(0.0);

  // AmanziMesh::Entity_kind kind = bc.kind();
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
    AmanziMesh::Entity_kind itkind;

    for (auto it = op->schema_row_.begin(); it != op->schema_row_.end(); ++it, ++item) {
      std::tie(itkind, std::ignore, std::ignore) = *it;

      if (itkind == AmanziMesh::NODE) {
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
      std::tie(itkind, std::ignore, std::ignore) = *it;

      if (itkind == AmanziMesh::NODE) {
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

              if (essential_eqn) {
                mesh_->node_get_cells(v, AmanziMesh::Parallel_type::ALL, &cells);
                rhs_loc(noff) = 0.0;
                if (v < nnodes_owned) (*rhs_node)[k][v] = value[k];
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


/* ******************************************************************
* Apply BCs of vector type.
****************************************************************** */
void PDE_HelperDiscretization::ApplyBCs_Cell_Vector_(
    const BCs& bc, Teuchos::RCP<Op> op,
    bool primary, bool eliminate, bool essential_eqn)
{
  const std::vector<int>& bc_model = bc.bc_model();
  const std::vector<std::vector<double> >& bc_value = bc.bc_value_vector();
  int d = bc_value[0].size();

  AmanziMesh::Entity_ID_List entities, cells;
  std::vector<int> offset;

  CompositeVector& rhs = *global_op_->rhs();
  rhs.PutScalarGhosted(0.0);

  const Schema& schema_row = global_op_->schema_row();
  const Schema& schema_col = global_op_->schema_col();

  AmanziMesh::Entity_kind kind = bc.kind();
  AMANZI_ASSERT(kind == AmanziMesh::FACE || kind == AmanziMesh::EDGE);
  Teuchos::RCP<Epetra_MultiVector> rhs_kind;
  if (primary) rhs_kind = rhs.ViewComponent(schema_row.KindToString(kind), true);

  for (int c = 0; c != ncells_owned; ++c) {
    WhetStone::DenseMatrix& Acell = op->matrices[c];
    int ncols = Acell.NumCols();
    int nrows = Acell.NumRows();

    int nents_owned(0);
    if (kind == AmanziMesh::FACE) {
      mesh_->cell_get_faces(c, &entities);
      nents_owned = ncells_owned;
    } else if (kind == AmanziMesh::EDGE) {
      mesh_->cell_get_edges(c, &entities);
      nents_owned = nedges_owned;
    }
    int nents = entities.size();

    // check for a boundary face
    bool found(false);
    for (int n = 0; n != nents; ++n) {
      int f = entities[n];
      if (bc_model[f] == OPERATOR_BC_DIRICHLET) found = true;
    }
    if (!found) continue;

    // essential conditions for test functions
    schema_row.ComputeOffset(c, mesh_, offset);

    bool flag(true);
    int item(0);
    AmanziMesh::Entity_kind itkind;

    for (auto it = op->schema_row_.begin(); it != op->schema_row_.end(); ++it, ++item) {
      std::tie(itkind, std::ignore, std::ignore) = *it;

      if (itkind == kind) {
        for (int n = 0; n != nents; ++n) {
          int x = entities[n];
          if (bc_model[x] == OPERATOR_BC_DIRICHLET) {
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
      std::tie(itkind, std::ignore, std::ignore) = *it;

      if (itkind == kind) {
        for (int n = 0; n != nents; ++n) {
          int x = entities[n];
          const std::vector<double>& value = bc_value[x];

          if (bc_model[x] == OPERATOR_BC_DIRICHLET) {
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

              if (essential_eqn) {
                rhs_loc(noff) = 0.0;
                (*rhs_kind)[k][x] = (x < nents_owned) ? value[k] : 0.0;

                if (kind == AmanziMesh::FACE) {
                  mesh_->face_get_cells(x, AmanziMesh::Parallel_type::ALL, &cells);
                } else if (kind == AmanziMesh::EDGE) {
                  mesh_->edge_get_cells(x, AmanziMesh::Parallel_type::ALL, &cells);
                }
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


/* ******************************************************************
* Composite vector space with one face component having multiple DOFs 
* and one regular cell component
****************************************************************** */
Teuchos::RCP<CompositeVectorSpace> CreateFracturedMatrixCVS(
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
    const Teuchos::RCP<const AmanziMesh::Mesh>& fracture)
{
  AmanziMesh::Entity_ID_List cells;
  int ncells_f = fracture->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  auto points = Teuchos::rcp(new Epetra_IntVector(mesh->face_map(true)));
  points->PutValue(1);

  for (int c = 0; c < ncells_f; ++c) {
    int f = fracture->entity_get_parent(AmanziMesh::CELL, c);
    mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    (*points)[f] = cells.size();
  }

  ParallelCommunication pp(mesh);
  pp.CopyMasterFace2GhostFace(*points);

  // create ghosted map with two points on each fracture face
  auto& gfmap = mesh->face_map(true);
  int nlocal = gfmap.NumMyElements();

  std::vector<int> gids(nlocal);
  gfmap.MyGlobalElements(&gids[0]);

  int* data; 
  points->ExtractView(&data);
  auto gmap = Teuchos::rcp(new Epetra_BlockMap(-1, nlocal, &gids[0], data, 0, gfmap.Comm()));

  // create master map with two points on each fracture face
  auto& mfmap = mesh->face_map(false);
  nlocal = mfmap.NumMyElements();

  auto mmap = Teuchos::rcp(new Epetra_BlockMap(-1, nlocal, &gids[0], data, 0, mfmap.Comm()));

  // create meta data
  std::string compname("face");
  auto cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh)->SetGhosted(true);
  cvs->AddComponent(compname, AmanziMesh::FACE, mmap, gmap, 1);
  cvs->AddComponent("cell", AmanziMesh::CELL, 1);

  return cvs;
}


/* ******************************************************************
* Composite vector space with one face component having multiple DOFs
****************************************************************** */
Teuchos::RCP<CompositeVectorSpace> CreateNonManifoldCVS(
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
{
  AmanziMesh::Entity_ID_List cells;
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);

  auto points = Teuchos::rcp(new Epetra_IntVector(mesh->face_map(true)));
  points->PutValue(1);

  for (int f = 0; f < nfaces; ++f) {
    mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    (*points)[f] = cells.size();
  }

  ParallelCommunication pp(mesh);
  pp.CopyMasterFace2GhostFace(*points);

  // create ghosted map with multiple points on each fracture face
  auto& gfmap = mesh->face_map(true);
  int nlocal = gfmap.NumMyElements();

  std::vector<int> gids(nlocal);
  gfmap.MyGlobalElements(&gids[0]);

  int* data; 
  points->ExtractView(&data);
  auto gmap = Teuchos::rcp(new Epetra_BlockMap(-1, nlocal, &gids[0], data, 0, gfmap.Comm()));

  // create master map with multiple points on each fracture face
  auto& mfmap = mesh->face_map(false);
  nlocal = mfmap.NumMyElements();

  auto mmap = Teuchos::rcp(new Epetra_BlockMap(-1, nlocal, &gids[0], data, 0, mfmap.Comm()));

  // create meta data
  std::string compname("face");
  auto cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh)->SetGhosted(true);
  cvs->AddComponent(compname, AmanziMesh::FACE, mmap, gmap, 1);

  return cvs;
}

}  // namespace Operators
}  // namespace Amanzi



