/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Helper class for discrete PDE operators.
*/

#include "errors.hh"

#include "Mesh_Algorithms.hh"
#include "WhetStoneDefs.hh"

#include "SchemaUtils.hh"
#include "ParallelCommunication.hh"
#include "PDE_HelperDiscretization.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Simple constructors.
****************************************************************** */
PDE_HelperDiscretization::PDE_HelperDiscretization(const Teuchos::RCP<Operator>& global_op)
  : global_op_(global_op)
{
  if (global_op == Teuchos::null) {
    Errors::Message msg("PDE_HelperDiscretization: Constructor received null global operator");
    Exceptions::amanzi_throw(msg);
  }

  mesh_ = global_op_->Mesh();
  PopulateDimensions_();
}


PDE_HelperDiscretization::PDE_HelperDiscretization(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  : mesh_(mesh)
{
  PopulateDimensions_();
}


PDE_HelperDiscretization::PDE_HelperDiscretization(const Teuchos::RCP<AmanziMesh::Mesh>& mesh)
  : mesh_(mesh)
{
  PopulateDimensions_();
}


/* ******************************************************************
* Supporting private routines.
****************************************************************** */
void
PDE_HelperDiscretization::PopulateDimensions_()
{
  ncells_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  nfaces_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  nnodes_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);

  ncells_wghost = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  nfaces_wghost = mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
  nnodes_wghost = mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL);

  if (mesh_->hasEdges()) {
    nedges_owned = mesh_->getNumEntities(AmanziMesh::Entity_kind::EDGE, AmanziMesh::Parallel_kind::OWNED);
    nedges_wghost = mesh_->getNumEntities(AmanziMesh::Entity_kind::EDGE, AmanziMesh::Parallel_kind::ALL);
  }
}


/* ******************************************************************
* Replace container of local matrices with another container.
****************************************************************** */
void
PDE_HelperDiscretization::set_local_op(const Teuchos::RCP<Op>& op)
{
  if (global_op_.get()) {
    if (local_op_.get()) {
      auto index =
        std::find(global_op_->begin(), global_op_->end(), local_op_) - global_op_->begin();
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
void
PDE_HelperDiscretization::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  auto base = global_op_->schema_row().get_base();

  for (auto bc : bcs_trial_) {
    bool missing(true);

    if (bc->type() == WhetStone::DOF_Type::SCALAR ||
        bc->type() == WhetStone::DOF_Type::NORMAL_COMPONENT ||
        bc->type() == WhetStone::DOF_Type::MOMENT) {
      if (base == AmanziMesh::Entity_kind::CELL) {
        ApplyBCs_Cell_Scalar_(*bc, local_op_, primary, eliminate, essential_eqn);
        missing = false;
      }
    } else if (bc->type() == WhetStone::DOF_Type::POINT) {
      if (base == AmanziMesh::Entity_kind::CELL) {
        ApplyBCs_Cell_Point_(*bc, local_op_, primary, eliminate, essential_eqn);
        missing = false;
      } else if (base == AmanziMesh::Entity_kind::NODE) {
        missing = false; // for testing
      }
    } else if (bc->type() == WhetStone::DOF_Type::VECTOR) {
      if (base == AmanziMesh::Entity_kind::CELL) {
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
void
PDE_HelperDiscretization::ApplyBCs_Cell_Scalar_(const BCs& bc,
                                                Teuchos::RCP<Op> op,
                                                bool primary,
                                                bool eliminate,
                                                bool essential_eqn)
{
  const std::vector<int>& bc_model = bc.bc_model();
  const std::vector<double>& bc_value = bc.bc_value();

  AmanziMesh::cEntity_ID_View entities, cells;
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
    if (kind == AmanziMesh::Entity_kind::FACE) {
      entities = mesh_->getCellFaces(c);
      nents_owned = nfaces_owned;
    } else if (kind == AmanziMesh::Entity_kind::EDGE) {
      entities = mesh_->getCellEdges(c);
      nents_owned = nedges_owned;
    } else if (kind == AmanziMesh::Entity_kind::NODE) {
      entities = mesh_->getCellNodes(c);
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

    for (auto it = op->schema_row().begin(); it != op->schema_row().end(); ++it, ++item) {
      std::tie(itkind, std::ignore, std::ignore) = *it;

      if (itkind == kind) {
        for (int n = 0; n != nents; ++n) {
          int x = entities[n];
          if (bc_model[x] == OPERATOR_BC_DIRICHLET) {
            if (flag) { // make a copy of elemental matrix
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
    for (auto it = op->schema_col().begin(); it != op->schema_col().end(); ++it, ++item) {
      std::tie(itkind, std::ignore, std::ignore) = *it;

      if (itkind == kind) {
        for (int n = 0; n != nents; ++n) {
          int x = entities[n];
          double value = bc_value[x];

          if (bc_model[x] == OPERATOR_BC_DIRICHLET) {
            if (flag) { // make a copy of elemental matrix
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

              if (kind == AmanziMesh::Entity_kind::FACE) {
                cells = mesh_->getFaceCells(x, AmanziMesh::Parallel_kind::ALL);
              } else if (kind == AmanziMesh::Entity_kind::NODE) {
                cells = mesh_->getNodeCells(x, AmanziMesh::Parallel_kind::ALL);
              } else if (kind == AmanziMesh::Entity_kind::EDGE) {
                cells = mesh_->getEdgeCells(x, AmanziMesh::Parallel_kind::ALL);
              }
              Acell(noff, noff) = 1.0 / cells.size();
            }

            AssembleVectorCellOp(c, *mesh_, schema_row, rhs_loc, rhs);
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
void
PDE_HelperDiscretization::ApplyBCs_Cell_Point_(const BCs& bc,
                                               Teuchos::RCP<Op> op,
                                               bool primary,
                                               bool eliminate,
                                               bool essential_eqn)
{
  const std::vector<int>& bc_model = bc.bc_model();
  const std::vector<AmanziGeometry::Point>& bc_value = bc.bc_value_point();
  std::vector<int> offset;

  CompositeVector& rhs = *global_op_->rhs();
  rhs.PutScalarGhosted(0.0);

  // AmanziMesh::Entity_kind kind = bc.kind();
  Teuchos::RCP<Epetra_MultiVector> rhs_node;
  if (primary) rhs_node = rhs.ViewComponent("node", true);

  int d = mesh_->getSpaceDimension();
  const Schema& schema_row = global_op_->schema_row();
  const Schema& schema_col = global_op_->schema_col();

  for (int c = 0; c != ncells_owned; ++c) {
    WhetStone::DenseMatrix& Acell = op->matrices[c];
    int ncols = Acell.NumCols();
    int nrows = Acell.NumRows();

    auto nodes = mesh_->getCellNodes(c);
    int nnodes = nodes.size();

    // check for a boundary face
    bool found(false);
    for (int n = 0; n != nnodes; ++n) {
      int v = nodes[n];
      if (bc_model[v] == OPERATOR_BC_DIRICHLET) found = true;
    }
    if (!found) continue;

    // essential conditions for test functions
    op->schema_row().ComputeOffset(c, mesh_, offset);

    bool flag(true);
    int item(0);
    AmanziMesh::Entity_kind itkind;

    for (auto it = op->schema_row().begin(); it != op->schema_row().end(); ++it, ++item) {
      std::tie(itkind, std::ignore, std::ignore) = *it;

      if (itkind == AmanziMesh::Entity_kind::NODE) {
        for (int n = 0; n != nnodes; ++n) {
          int v = nodes[n];
          if (bc_model[v] == OPERATOR_BC_DIRICHLET) {
            if (flag) { // make a copy of elemental matrix
              op->matrices_shadow[c] = Acell;
              flag = false;
            }
            for (int k = 0; k < d; ++k) {
              int noff(d * n + k + offset[item]);
              for (int m = 0; m < ncols; m++) Acell(noff, m) = 0.0;
            }
          }
        }
      }
    }

    // essential zero conditions for trial functions
    schema_col.ComputeOffset(c, mesh_, offset);

    item = 0;
    for (auto it = op->schema_col().begin(); it != op->schema_col().end(); ++it, ++item) {
      std::tie(itkind, std::ignore, std::ignore) = *it;

      if (itkind == AmanziMesh::Entity_kind::NODE) {
        for (int n = 0; n != nnodes; ++n) {
          int v = nodes[n];
          AmanziGeometry::Point value = bc_value[v];

          if (bc_model[v] == OPERATOR_BC_DIRICHLET) {
            if (flag) { // make a copy of elemental matrix
              op->matrices_shadow[c] = Acell;
              flag = false;
            }

            for (int k = 0; k < d; ++k) {
              int noff(d * n + k + offset[item]);
              WhetStone::DenseVector rhs_loc(nrows);

              if (eliminate) {
                for (int m = 0; m < nrows; m++) {
                  rhs_loc(m) = -Acell(m, noff) * value[k];
                  Acell(m, noff) = 0.0;
                }
              }

              if (essential_eqn) {
                auto cells = mesh_->getNodeCells(v, AmanziMesh::Parallel_kind::ALL);
                rhs_loc(noff) = 0.0;
                if (v < nnodes_owned) (*rhs_node)[k][v] = value[k];
                Acell(noff, noff) = 1.0 / cells.size();
              }

              AssembleVectorCellOp(c, *mesh_, schema_row, rhs_loc, rhs);
            }
          }
        }
      }
    }
  }

  rhs.GatherGhostedToMaster(Add);
}


/* ******************************************************************
* Apply BCs of vector type. The code is based on face (f) DOFs.
****************************************************************** */
void
PDE_HelperDiscretization::ApplyBCs_Cell_Vector_(const BCs& bc,
                                                Teuchos::RCP<Op> op,
                                                bool primary,
                                                bool eliminate,
                                                bool essential_eqn)
{
  const std::vector<int>& bc_model = bc.bc_model();
  const std::vector<std::vector<double>>& bc_value = bc.bc_value_vector();
  int d = bc_value[0].size();

  AmanziMesh::cEntity_ID_View entities;
  std::vector<int> offset;

  CompositeVector& rhs = *global_op_->rhs();
  rhs.PutScalarGhosted(0.0);

  const Schema& schema_row = global_op_->schema_row();
  const Schema& schema_col = global_op_->schema_col();

  AmanziMesh::Entity_kind kind = bc.kind();
  AMANZI_ASSERT(kind == AmanziMesh::Entity_kind::FACE || kind == AmanziMesh::Entity_kind::EDGE);
  Teuchos::RCP<Epetra_MultiVector> rhs_kind;
  if (primary) rhs_kind = rhs.ViewComponent(schema_row.KindToString(kind), true);

  for (int c = 0; c != ncells_owned; ++c) {
    WhetStone::DenseMatrix& Acell = op->matrices[c];
    int ncols = Acell.NumCols();
    int nrows = Acell.NumRows();

    if (kind == AmanziMesh::Entity_kind::FACE) { entities = mesh_->getCellFaces(c); }
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

    for (auto it = op->schema_row().begin(); it != op->schema_row().end(); ++it, ++item) {
      std::tie(itkind, std::ignore, std::ignore) = *it;

      if (itkind == kind) {
        for (int n = 0; n != nents; ++n) {
          int f = entities[n];
          if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
            if (flag) { // make a copy of elemental matrix
              op->matrices_shadow[c] = Acell;
              flag = false;
            }

            for (int k = 0; k < d; ++k) {
              int noff(d * n + k + offset[item]);
              for (int m = 0; m < ncols; m++) Acell(noff, m) = 0.0;
            }
          }
        }
      }
    }

    // essential zero conditions for trial functions
    schema_col.ComputeOffset(c, mesh_, offset);

    item = 0;
    for (auto it = op->schema_col().begin(); it != op->schema_col().end(); ++it, ++item) {
      std::tie(itkind, std::ignore, std::ignore) = *it;

      if (itkind == kind) {
        for (int n = 0; n != nents; ++n) {
          int f = entities[n];
          const std::vector<double>& value = bc_value[f];

          if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
            if (flag) { // make a copy of elemental matrix
              op->matrices_shadow[c] = Acell;
              flag = false;
            }

            for (int k = 0; k < d; ++k) {
              int noff(d * n + k + offset[item]);
              WhetStone::DenseVector rhs_loc(nrows);

              if (eliminate) {
                for (int m = 0; m < nrows; m++) {
                  rhs_loc(m) = -Acell(m, noff) * value[k];
                  Acell(m, noff) = 0.0;
                }
              }

              if (essential_eqn) {
                rhs_loc(noff) = 0.0;
                (*rhs_kind)[k][f] = value[k];
                Acell(noff, noff) = 1.0;
              }

              AssembleVectorCellOp(c, *mesh_, schema_row, rhs_loc, rhs);
            }
          }
        }
      }
    }
  }

  rhs.GatherGhostedToMaster(Add);
}


/* ******************************************************************
* Enforce essential boundary conditions.
****************************************************************** */
void
PDE_HelperDiscretization::EnforceBCs(CompositeVector& field)
{
  for (auto bc : bcs_trial_) {
    std::string name = to_string(bc->kind());
    if (field.HasComponent(name)) {
      auto& field_comp = *field.ViewComponent(name, true);

      const std::vector<int>& bc_model = bc->bc_model();

      if (bc->type() == WhetStone::DOF_Type::SCALAR ||
          bc->type() == WhetStone::DOF_Type::NORMAL_COMPONENT) {
        const auto& bc_value = bc->bc_value();
        for (int i = 0; i < bc_model.size(); ++i) {
          if (bc_model[i] == OPERATOR_BC_DIRICHLET) { field_comp[0][i] = bc_value[i]; }
        }
      } else if (bc->type() == WhetStone::DOF_Type::VECTOR) {
        const auto& bc_value = bc->bc_value_vector();
        int n = bc_value[0].size();
        for (int i = 0; i < bc_model.size(); ++i) {
          if (bc_model[i] == OPERATOR_BC_DIRICHLET) {
            for (int k = 0; k < n; ++k) field_comp[k][i] = bc_value[i][k];
          }
        }
      } else if (bc->type() == WhetStone::DOF_Type::POINT) {
        const auto& bc_value = bc->bc_value_point();
        int n = bc_value[0].dim();
        for (int i = 0; i < bc_model.size(); ++i) {
          if (bc_model[i] == OPERATOR_BC_DIRICHLET) {
            for (int k = 0; k < n; ++k) field_comp[k][i] = bc_value[i][k];
          }
        }
      }
    }
  }
}


/* ******************************************************************
* Composite vector space with one face component having multiple DOFs
* and one regular cell component
****************************************************************** */
Teuchos::RCP<CompositeVectorSpace>
CreateFracturedMatrixCVS(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                         const Teuchos::RCP<const AmanziMesh::Mesh>& fracture)
{
  int ncells_f = fracture->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  auto points = Teuchos::rcp(new Epetra_IntVector(mesh->getMap(AmanziMesh::Entity_kind::FACE,true)));
  points->PutValue(1);

  for (int c = 0; c < ncells_f; ++c) {
    int f = fracture->getEntityParent(AmanziMesh::Entity_kind::CELL, c);
    auto cells = mesh->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    (*points)[f] = cells.size();
  }

  ParallelCommunication pp(mesh);
  pp.CopyMasterFace2GhostFace(*points);

  // create ghosted map with two points on each fracture face
  auto& gfmap = mesh->getMap(AmanziMesh::Entity_kind::FACE,true);
  int nlocal = gfmap.NumMyElements();

  std::vector<int> gids(nlocal);
  gfmap.MyGlobalElements(&gids[0]);

  int* data;
  points->ExtractView(&data);
  auto gmap = Teuchos::rcp(new Epetra_BlockMap(-1, nlocal, &gids[0], data, 0, gfmap.Comm()));

  // create master map with two points on each fracture face
  auto& mfmap = mesh->getMap(AmanziMesh::Entity_kind::FACE,false);
  nlocal = mfmap.NumMyElements();

  auto mmap = Teuchos::rcp(new Epetra_BlockMap(-1, nlocal, &gids[0], data, 0, mfmap.Comm()));

  // create meta data
  std::string compname("face");
  auto cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh)->SetGhosted(true);
  cvs->AddComponent(compname, AmanziMesh::Entity_kind::FACE, mmap, gmap, 1);
  cvs->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  return cvs;
}


/* ******************************************************************
* Composite vector space with some faces having multiple DOFs
****************************************************************** */
Teuchos::RCP<CompositeVectorSpace>
CreateManifoldCVS(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
{
  int nfaces = mesh->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);

  auto points = Teuchos::rcp(new Epetra_IntVector(mesh->getMap(AmanziMesh::Entity_kind::FACE,true)));
  points->PutValue(1);

  for (int f = 0; f < nfaces; ++f) {
   auto  cells = mesh->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    (*points)[f] = cells.size();
  }

  ParallelCommunication pp(mesh);
  pp.CopyMasterFace2GhostFace(*points);

  // create ghosted map with multiple points on each fracture face
  auto& gfmap = mesh->getMap(AmanziMesh::Entity_kind::FACE,true);
  int nlocal = gfmap.NumMyElements();

  std::vector<int> gids(nlocal);
  gfmap.MyGlobalElements(&gids[0]);

  int* data;
  points->ExtractView(&data);
  auto gmap = Teuchos::rcp(new Epetra_BlockMap(-1, nlocal, &gids[0], data, 0, gfmap.Comm()));

  // create master map with multiple points on each fracture face
  auto& mfmap = mesh->getMap(AmanziMesh::Entity_kind::FACE,false);
  nlocal = mfmap.NumMyElements();

  auto mmap = Teuchos::rcp(new Epetra_BlockMap(-1, nlocal, &gids[0], data, 0, mfmap.Comm()));

  // create meta data
  std::string compname("face");
  auto cvs = Teuchos::rcp(new CompositeVectorSpace());
  cvs->SetMesh(mesh)->SetGhosted(true);
  cvs->AddComponent(compname, AmanziMesh::Entity_kind::FACE, mmap, gmap, 1);

  return cvs;
}


/* ******************************************************************
* Support function: copy data from cells to dirichlet faces
****************************************************************** */
void
CellToBoundaryFaces(const std::vector<int>& bc_model, CompositeVector& field)
{
  int nfaces = bc_model.size();
  const auto& mesh = field.Mesh();
  auto& field_c = *field.ViewComponent("cell", true);

  if (field.HasComponent("face")) {
    const auto& fmap = *field.ComponentMap("face", true);
    auto& field_f = *field.ViewComponent("face", true);

    for (int f = 0; f != nfaces; ++f) {
      if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
        int g = fmap.FirstPointInElement(f);
        int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh, f);
        field_f[0][g] = field_c[0][c];
      }
    }
  } else {
    auto& field_bf = *field.ViewComponent("boundary_face", true);

    for (int f = 0; f != nfaces; ++f) {
      if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
        int bf = AmanziMesh::getFaceOnBoundaryBoundaryFace(*mesh, f);
        int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh, f);
        field_bf[0][bf] = field_c[0][c];
      }
    }
  }
}


/* ******************************************************************
* Support function: copy data from boundary faces to faces (if any)
****************************************************************** */
void
BoundaryFacesToFaces(const std::vector<int>& bc_model,
                     const CompositeVector& input,
                     CompositeVector& output)
{
  int nfaces = bc_model.size();
  const auto& mesh = output.Mesh();

  const auto& fmap = *output.ComponentMap("face", true);
  const auto& input_bf = *input.ViewComponent("boundary_face", true);
  auto& output_f = *output.ViewComponent("face", true);

  for (int f = 0; f != nfaces; ++f) {
    if (bc_model[f] != Operators::OPERATOR_BC_NONE) {
      int g = fmap.FirstPointInElement(f);
      int bf = AmanziMesh::getFaceOnBoundaryBoundaryFace(*mesh, f);
      output_f[0][g] = input_bf[0][bf];
    }
  }
}


/* ******************************************************************
* Support function: copy data from BCs to faces
****************************************************************** */
void
BoundaryDataToFaces(const Teuchos::RCP<Operators::BCs>& op_bc, CompositeVector& field)
{
  std::vector<int>& bc_model = op_bc->bc_model();
  std::vector<double>& bc_value = op_bc->bc_value();

  int nfaces = bc_model.size();
  const auto& mesh = field.Mesh();

  if (field.HasComponent("face")) {
    const auto& fmap = *field.ComponentMap("face", true);
    auto& field_f = *field.ViewComponent("face", true);

    for (int f = 0; f != nfaces; ++f) {
      if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) { 
        int g = fmap.FirstPointInElement(f);
        field_f[0][g] = bc_value[f];
      }
    }
  } else {
    auto& field_bf = *field.ViewComponent("boundary_face", true);

    for (int f = 0; f != nfaces; ++f) {
      if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
        int bf = AmanziMesh::getFaceOnBoundaryBoundaryFace(*mesh, f);
        field_bf[0][bf] = bc_value[f];
      }
    }
  }
}

} // namespace Operators
} // namespace Amanzi
