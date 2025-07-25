/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)
*/

/*
  Operators

  This operator is a collection of local "DIAGONAL" Ops.
*/

#include "WhetStoneMeshUtils.hh"

#include "OperatorUtils.hh"
#include "Operator_Cell.hh"
#include "Operator_Edge.hh"
#include "Operator_Node.hh"
#include "Op_Cell_Cell.hh"
#include "Op_Face_Face.hh"
#include "Op_Edge_Edge.hh"
#include "Op_Node_Node.hh"
#include "Op_SurfaceCell_SurfaceCell.hh"
#include "PDE_Accumulation.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Modifier for diagonal operators.  Op += du
****************************************************************** */
void
PDE_Accumulation::AddAccumulationTerm(const CompositeVector& du, const std::string& name)
{
  Teuchos::RCP<Op> op = FindOp_(name);
  Epetra_MultiVector& diag = *op->diag;
  const Epetra_MultiVector& duc = *du.ViewComponent(name);

  int n = duc.MyLength();
  int m = duc.NumVectors();
  int m0 = diag.NumVectors();

  if (m == m0) {
    for (int k = 0; k < m; k++) {
      for (int i = 0; i < n; i++) {
        diag[k][i] += duc[k][i];
      }
    }
  } else if (m == 1) {
    for (int k = 0; k < m0; k++) {
      for (int i = 0; i < n; i++) {
        diag[k][i] += duc[0][i];
      }
    }
  } else {
    AMANZI_ASSERT(false);
  }
}


/* ******************************************************************
* Modifier for diagonal operators.  Op += du * vol / dt.
****************************************************************** */
void
PDE_Accumulation::AddAccumulationTerm(const CompositeVector& du,
                                      double dT,
                                      const std::string& name,
                                      bool volume)
{
  Teuchos::RCP<Op> op = FindOp_(name);
  Epetra_MultiVector& diag = *op->diag;

  const Epetra_MultiVector& duc = *du.ViewComponent(name);

  int n = duc.MyLength();
  int m = duc.NumVectors();

  if (volume) {
    CompositeVector vol(du);
    CalculateEntityVolume_(vol, name);
    Epetra_MultiVector& volc = *vol.ViewComponent(name);

    for (int k = 0; k < m; k++) {
      for (int i = 0; i < n; i++) diag[k][i] += volc[0][i] * duc[k][i] / dT;
    }

  } else {
    for (int k = 0; k < m; k++) {
      for (int i = 0; i < n; i++) diag[k][i] += duc[k][i] / dT;
    }
  }
}


/* ******************************************************************
* Op += alpha * du1 * du2 * vol. Both du1 and du2 may scalar.
****************************************************************** */
void
PDE_Accumulation::AddAccumulationTerm(const CompositeVector& du1,
                                      const CompositeVector& du2,
                                      double alpha,
                                      const std::string& name,
                                      bool volume)
{
  Teuchos::RCP<Op> op = FindOp_(name);
  Epetra_MultiVector& diag = *op->diag;

  const Epetra_MultiVector& du1c = *du1.ViewComponent(name);
  const Epetra_MultiVector& du2c = *du2.ViewComponent(name);

  int n1 = du1c.MyLength();
  int m1 = du1c.NumVectors();

  int n2 = du2c.MyLength();
  int m2 = du2c.NumVectors();

  int n0 = diag.MyLength();
  int m0 = diag.NumVectors();

  AMANZI_ASSERT(m1 == m0 || m1 == 1);
  AMANZI_ASSERT(m2 == m0 || m2 == 1);
  AMANZI_ASSERT(n1 == n0 && n2 == n0);

  if (volume) {
    CompositeVector vol(du1);
    CalculateEntityVolume_(vol, name);
    Epetra_MultiVector& volc = *vol.ViewComponent(name);

    for (int k = 0; k < m0; k++) {
      int k0 = std::min(0, k);
      for (int i = 0; i < n0; i++) diag[k][i] += alpha * volc[0][i] * du1c[k0][i] * du2c[k0][i];
    }

  } else {
    for (int k = 0; k < m0; k++) {
      int k0 = std::min(0, k);
      for (int i = 0; i < n0; i++) {
        diag[k][i] += alpha * du1c[k0][i] * du2c[k0][i];
      }
    }
  }
}


/* ******************************************************************
* Modifier for diagonal operators and rhs.
* Op  += alpha * s1 * vol
* Rhs += alpha * s2 * vol
****************************************************************** */
void
PDE_Accumulation::AddAccumulationRhs(const CompositeVector& s1,
                                     const CompositeVector& s2,
                                     double alpha,
                                     const std::string& name,
                                     bool volume)
{
  Teuchos::RCP<Op> op = FindOp_(name);
  Epetra_MultiVector& diag = *op->diag;

  const Epetra_MultiVector& s1c = *s1.ViewComponent(name);
  const Epetra_MultiVector& s2c = *s2.ViewComponent(name);

  int n = s1c.MyLength();
  int m = s1c.NumVectors();

  Epetra_MultiVector& rhs = *global_operator()->rhs()->ViewComponent(name);

  AMANZI_ASSERT(s1c.MyLength() == s2c.MyLength());
  AMANZI_ASSERT(s1c.MyLength() == diag.MyLength());
  AMANZI_ASSERT(s2c.MyLength() == rhs.MyLength());

  AMANZI_ASSERT(s1c.NumVectors() == s2c.NumVectors());
  AMANZI_ASSERT(s1c.NumVectors() == diag.NumVectors());
  AMANZI_ASSERT(s2c.NumVectors() == rhs.NumVectors());

  if (volume) {
    CompositeVector vol(s1);
    CalculateEntityVolume_(vol, name);
    Epetra_MultiVector& volc = *vol.ViewComponent(name);

    for (int k = 0; k < m; k++) {
      for (int i = 0; i < n; i++) {
        diag[k][i] += volc[0][i] * s1c[k][i] * alpha;
        rhs[k][i] += volc[0][i] * s2c[k][i] * alpha;
      }
    }

  } else {
    for (int k = 0; k < m; k++) {
      for (int i = 0; i < n; i++) {
        diag[k][i] += s1c[k][i] * alpha;
        rhs[k][i] += s2c[k][i] * alpha;
      }
    }
  }
}


/* ******************************************************************
* Linearized update methods with storage terms for component "name".
* Op  += ss * vol / dt
* RHS += s0 * vol * u0 / dt
****************************************************************** */
void
PDE_Accumulation::AddAccumulationDelta(const CompositeVector& u0,
                                       const CompositeVector& s0,
                                       const CompositeVector& ss,
                                       double dT,
                                       const std::string& name)
{
  Teuchos::RCP<Op> op = FindOp_(name);
  Epetra_MultiVector& diag = *op->diag;

  CompositeVector vol(ss);
  CalculateEntityVolume_(vol, name);

  const Epetra_MultiVector& u0c = *u0.ViewComponent(name);
  const Epetra_MultiVector& s0c = *s0.ViewComponent(name);
  const Epetra_MultiVector& ssc = *ss.ViewComponent(name);

  Epetra_MultiVector& volc = *vol.ViewComponent(name);
  Epetra_MultiVector& rhs = *global_operator()->rhs()->ViewComponent(name);

  int n = u0c.MyLength();
  int m = u0c.NumVectors();
  for (int k = 0; k < m; ++k) {
    for (int i = 0; i < n; i++) {
      double factor = volc[0][i] / dT;
      diag[k][i] += factor * ssc[k][i];
      rhs[k][i] += factor * s0c[k][i] * u0c[k][i];
    }
  }
}


/* ******************************************************************
* Linearized update methods with storage terms for component "name".
* Op  += vol / dt
* RHS += vol * u0 / dt
****************************************************************** */
void
PDE_Accumulation::AddAccumulationDelta(const CompositeVector& u0,
                                       double dT,
                                       const std::string& name)
{
  Teuchos::RCP<Op> op = FindOp_(name);
  Epetra_MultiVector& diag = *op->diag;

  CompositeVector vol(u0);
  CalculateEntityVolume_(vol, name);

  const Epetra_MultiVector& u0c = *u0.ViewComponent(name);
  Epetra_MultiVector& volc = *vol.ViewComponent(name);
  Epetra_MultiVector& rhs = *global_operator()->rhs()->ViewComponent(name);

  int n = u0c.MyLength();
  int m = u0c.NumVectors();
  for (int k = 0; k < m; ++k) {
    for (int i = 0; i < n; i++) {
      double factor = volc[0][i] / dT;
      diag[k][i] += factor;
      rhs[k][i] += factor * u0c[k][i];
    }
  }
}


/* ******************************************************************
* Linearized update methods with storage terms for component "name".
* Op  += ss
* RHS += ss * u0
****************************************************************** */
void
PDE_Accumulation::AddAccumulationDeltaNoVolume(const CompositeVector& u0,
                                               const CompositeVector& ss,
                                               const std::string& name)
{
  if (!ss.HasComponent(name) ) AMANZI_ASSERT(false);

  Teuchos::RCP<Op> op = FindOp_(name);
  Epetra_MultiVector& diag = *op->diag;

  const Epetra_MultiVector& u0c = *u0.ViewComponent(name);
  const Epetra_MultiVector& ssc = *ss.ViewComponent(name);

  Epetra_MultiVector& rhs = *global_operator()->rhs()->ViewComponent(name);

  int n = u0c.MyLength();
  int m = u0c.NumVectors();
  for (int k = 0; k < m; ++k) {
    for (int i = 0; i < n; i++) {
      diag[k][i] += ssc[k][i];
      rhs[k][i] += ssc[k][i] * u0c[k][i];
    }
  }
}


/* ******************************************************************
* Calculate entity volume for component "name" of vector ss.
****************************************************************** */
void
PDE_Accumulation::CalculateEntityVolume_(CompositeVector& volume, const std::string& name)
{
  if (name == "cell" && volume.HasComponent("cell")) {
    Epetra_MultiVector& vol = *volume.ViewComponent(name);

    for (int c = 0; c != ncells_owned; ++c) {
      vol[0][c] = mesh_->getCellVolume(c);
    }

  } else if (name == "face" && volume.HasComponent("face")) {
    // Missing code.
    AMANZI_ASSERT(false);

  } else if (name == "edge" && volume.HasComponent("edge")) {
    Epetra_MultiVector& vol = *volume.ViewComponent(name, true);
    vol.PutScalar(0.0);

    for (int c = 0; c != ncells_owned; ++c) {
      auto edges = mesh_->getCellEdges(c);
      int nedges = edges.size();

      for (int i = 0; i < nedges; i++) {
        vol[0][edges[i]] += mesh_->getCellVolume(c) / nedges;
      }
    }
    volume.GatherGhostedToMaster(name);

  } else if (name == "node" && volume.HasComponent("node")) {
    Epetra_MultiVector& vol = *volume.ViewComponent(name, true);
    vol.PutScalar(0.0);

    for (int c = 0; c != ncells_owned; ++c) {
      auto nodes = mesh_->getCellNodes(c);
      int nnodes = nodes.size();

      double cellvolume = mesh_->getCellVolume(c);
      std::vector<double> weights(nnodes, 1.0 / nnodes);

      if (mesh_->getSpaceDimension() == 2) {
        WhetStone::PolygonCentroidWeights(*mesh_, nodes, cellvolume, weights);
      }

      for (int i = 0; i < nnodes; i++) {
        vol[0][nodes[i]] += weights[i] * cellvolume;
      }
    }
    volume.GatherGhostedToMaster(name);

  } else {
    AMANZI_ASSERT(false);
  }
}


/* ******************************************************************
* Note: When complex schema is used to create a set of local ops, the
* the local local_op_ is not well defined.
****************************************************************** */
void
PDE_Accumulation::Init_(const Schema& schema, bool surf)
{
  int num;
  AmanziMesh::Entity_kind kind;
  Errors::Message msg;

  if (global_op_ == Teuchos::null) {
    // constructor was given a mesh
    global_op_schema_ = schema;
    local_op_schema_ = schema;

    for (auto it = schema.begin(); it != schema.end(); ++it) {
      std::tie(kind, std::ignore, num) = *it;

      Teuchos::RCP<Op> op;
      auto cvs = CreateCompositeVectorSpace(mesh_, AmanziMesh::to_string(kind), kind, num);

      if (kind == AmanziMesh::Entity_kind::CELL) {
        int old_schema = OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_CELL;
        global_op_ = Teuchos::rcp(new Operator_Cell(cvs, plist_, old_schema));
        std::string name("CELL_CELL");
        if (surf) {
          op = Teuchos::rcp(new Op_SurfaceCell_SurfaceCell(name, mesh_));
        } else {
          op = Teuchos::rcp(new Op_Cell_Cell(name, mesh_, num));
        }

      } else if (kind == AmanziMesh::Entity_kind::EDGE) {
        global_op_ = Teuchos::rcp(new Operator_Edge(cvs, plist_));
        std::string name("EDGE_EDGE");
        op = Teuchos::rcp(new Op_Edge_Edge(name, mesh_));

      } else if (kind == AmanziMesh::Entity_kind::NODE) {
        global_op_ = Teuchos::rcp(new Operator_Node(cvs, plist_));
        std::string name("NODE_NODE");
        op = Teuchos::rcp(new Op_Node_Node(name, mesh_, num));

      } else {
        msg << "Accumulation operator: Unknown kind \"" << AmanziMesh::to_string(kind) << "\".\n";
        Exceptions::amanzi_throw(msg);
      }

      global_op_->OpPushBack(op);
      local_ops_.push_back(op);
    }

  } else {
    // constructor was given an Operator
    global_op_schema_ = global_op_->schema_row();
    mesh_ = global_op_->DomainMap().Mesh();

    for (auto it = schema.begin(); it != schema.end(); ++it) {
      std::tie(kind, std::ignore, num) = *it;

      int old_schema(0);
      Teuchos::RCP<Op> op;

      if (kind == AmanziMesh::Entity_kind::CELL) {
        old_schema = OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_CELL;
        std::string name("CELL_CELL");
        if (surf) {
          op = Teuchos::rcp(new Op_SurfaceCell_SurfaceCell(name, mesh_));
        } else {
          op = Teuchos::rcp(new Op_Cell_Cell(name, mesh_, num));
        }

      } else if (kind == AmanziMesh::Entity_kind::FACE) {
        old_schema = OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_FACE;
        std::string name("FACE_FACE");
        op = Teuchos::rcp(new Op_Face_Face(name, mesh_));

      } else if (kind == AmanziMesh::Entity_kind::EDGE) {
        old_schema = OPERATOR_SCHEMA_BASE_EDGE | OPERATOR_SCHEMA_DOFS_EDGE;
        std::string name("EDGE_EDGE");
        op = Teuchos::rcp(new Op_Edge_Edge(name, mesh_));

      } else if (kind == AmanziMesh::Entity_kind::NODE) {
        old_schema = OPERATOR_SCHEMA_BASE_NODE | OPERATOR_SCHEMA_DOFS_NODE;
        std::string name("NODE_NODE");
        op = Teuchos::rcp(new Op_Node_Node(name, mesh_, num));

      } else {
        msg << "Accumulation operator: Unknown kind \"" << AmanziMesh::to_string(kind) << "\".\n";
        Exceptions::amanzi_throw(msg);
      }

      // register the accumulation Op
      local_op_schema_.Init(old_schema);
      global_op_->OpPushBack(op);
      local_ops_.push_back(op);
    }
  }

  // mesh info
  ncells_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  nfaces_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  nnodes_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
}


/* ******************************************************************
* Apply boundary conditions to
****************************************************************** */
void
PDE_Accumulation::ApplyBCs()
{
  for (auto bc = bcs_trial_.begin(); bc != bcs_trial_.end(); ++bc) {
    const std::vector<int>& bc_model = (*bc)->bc_model();

    for (auto it = local_ops_.begin(); it != local_ops_.end(); ++it) {
      const Schema& schema = (*it)->schema_row();
      if (schema.get_base() == (*bc)->kind()) {
        Epetra_MultiVector& diag = *(*it)->diag;
        int m = diag.NumVectors();

        for (int i = 0; i < diag.MyLength(); i++) {
          if (bc_model[i] == OPERATOR_BC_DIRICHLET) {
            for (int k = 0; k < m; ++k) {
              diag[k][i] = 0.0;
            }
          } else if (bc_model[i] == OPERATOR_BC_KINEMATIC) {
            const auto& normal = WhetStone::getNodeUnitNormal(*mesh_, i);
            int k = (std::fabs(normal[0]) > std::fabs(normal[1])) ? 0 : 1;
            if (normal.dim() == 3) k = (std::fabs(normal[k]) > std::fabs(normal[2])) ? k : 2;

            diag[k][i] = 0.0;
          }
        }
      }
    }
  }
}


/* ******************************************************************
* Return operator with the given base.
****************************************************************** */
Teuchos::RCP<Op>
PDE_Accumulation::FindOp_(const std::string& name) const
{
  for (auto it = local_ops_.begin(); it != local_ops_.end(); ++it) {
    const Schema& schema = (*it)->schema_row();
    if (AmanziMesh::to_string(schema.get_base() ) == name) return *it;
  }
  return Teuchos::null;
}

} // namespace Operators
} // namespace Amanzi
