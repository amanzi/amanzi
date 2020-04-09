/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "WhetStoneMeshUtils.hh"

#include "OperatorUtils.hh"
#include "Operator_Cell.hh"
#include "Operator_Edge.hh"
#include "Operator_Node.hh"
#include "Op_Cell_Cell.hh"
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
PDE_Accumulation::AddAccumulationTerm(const CompositeVector& du,
                                      const std::string& name)
{
  Teuchos::RCP<Op> op = FindOp_(name);
  Epetra_MultiVector& diag = *op->diag;
  const Epetra_MultiVector& duc = *du.ViewComponent(name);

  int n = duc.getLocalLength();
  int m = duc.getNumVectors();
  for (int k = 0; k < m; k++) {
    for (int i = 0; i < n; i++) { diag[k][i] += duc[k][i]; }
  }
}


/* ******************************************************************
 * Modifier for diagonal operators.  Op += du * vol / dt.
 ****************************************************************** */
void
PDE_Accumulation::AddAccumulationTerm(const CompositeVector& du, double dT,
                                      const std::string& name, bool volume)
{
  Teuchos::RCP<Op> op = FindOp_(name);
  Epetra_MultiVector& diag = *op->diag;

  const Epetra_MultiVector& duc = *du.ViewComponent(name);

  int n = duc.getLocalLength();
  int m = duc.getNumVectors();

  if (volume) {
    CompositeVector vol(du);
    CalculateEntityVolume_(vol, name);
    Epetra_MultiVector& volc = *vol.ViewComponent(name);

    for (int k = 0; k < m; k++) {
      for (int i = 0; i < n; i++) { diag[k][i] += volc[0][i] * duc[k][i] / dT; }
    }

  } else {
    for (int k = 0; k < m; k++) {
      for (int i = 0; i < n; i++) { diag[k][i] += duc[k][i] / dT; }
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
                                     const CompositeVector& s2, double alpha,
                                     const std::string& name, bool volume)
{
  Teuchos::RCP<Op> op = FindOp_(name);
  Epetra_MultiVector& diag = *op->diag;

  const Epetra_MultiVector& s1c = *s1.ViewComponent(name);
  const Epetra_MultiVector& s2c = *s2.ViewComponent(name);

  int n = s1c.getLocalLength();
  int m = s1c.getNumVectors();

  Epetra_MultiVector& rhs = *global_operator()->rhs()->ViewComponent(name);

  AMANZI_ASSERT(s1c.getLocalLength() == s2c.getLocalLength());
  AMANZI_ASSERT(s1c.getLocalLength() == diag.getLocalLength());
  AMANZI_ASSERT(s2c.getLocalLength() == rhs.getLocalLength());

  AMANZI_ASSERT(s1c.getNumVectors() == s2c.getNumVectors());
  AMANZI_ASSERT(s1c.getNumVectors() == diag.getNumVectors());
  AMANZI_ASSERT(s2c.getNumVectors() == rhs.getNumVectors());

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
                                       const CompositeVector& ss, double dT,
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

  int n = u0c.getLocalLength();
  int m = u0c.getNumVectors();
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
PDE_Accumulation::AddAccumulationDelta(const CompositeVector& u0, double dT,
                                       const std::string& name)
{
  Teuchos::RCP<Op> op = FindOp_(name);
  Epetra_MultiVector& diag = *op->diag;

  CompositeVector vol(u0);
  CalculateEntityVolume_(vol, name);

  const Epetra_MultiVector& u0c = *u0.ViewComponent(name);
  Epetra_MultiVector& volc = *vol.ViewComponent(name);
  Epetra_MultiVector& rhs = *global_operator()->rhs()->ViewComponent(name);

  int n = u0c.getLocalLength();
  int m = u0c.getNumVectors();
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
  if (!ss.HasComponent(name)) AMANZI_ASSERT(false);

  Teuchos::RCP<Op> op = FindOp_(name);
  Epetra_MultiVector& diag = *op->diag;

  const Epetra_MultiVector& u0c = *u0.ViewComponent(name);
  const Epetra_MultiVector& ssc = *ss.ViewComponent(name);

  Epetra_MultiVector& rhs = *global_operator()->rhs()->ViewComponent(name);

  int n = u0c.getLocalLength();
  int m = u0c.getNumVectors();
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
PDE_Accumulation::CalculateEntityVolume_(CompositeVector& volume,
                                         const std::string& name)
{
  AmanziMesh::Entity_ID_List nodes, edges;

  if (name == "cell" && volume.HasComponent("cell")) {
    Epetra_MultiVector& vol = *volume.ViewComponent(name);

    for (int c = 0; c != ncells_owned; ++c) {
      vol[0][c] = mesh_->cell_volume(c, false);
    }

  } else if (name == "face" && volume.HasComponent("face")) {
    // Missing code.
    AMANZI_ASSERT(false);

  } else if (name == "edge" && volume.HasComponent("edge")) {
    Epetra_MultiVector& vol = *volume.ViewComponent(name, true);
    vol.putScalar(0.0);

    for (int c = 0; c != ncells_owned; ++c) {
      mesh_->cell_get_edges(c, &edges);
      int nedges = edges.size();

      for (int i = 0; i < nedges; i++) {
        vol[0][edges[i]] += mesh_->cell_volume(c, false) / nedges;
      }
    }
    volume.GatherGhostedToMaster(name);

  } else if (name == "node" && volume.HasComponent("node")) {
    Epetra_MultiVector& vol = *volume.ViewComponent(name, true);
    vol.putScalar(0.0);

    for (int c = 0; c != ncells_owned; ++c) {
      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();

      double cellvolume = mesh_->cell_volume(c, false);
      std::vector<double> weights(nnodes, 1.0 / nnodes);

      if (mesh_->space_dimension() == 2) {
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
PDE_Accumulation::InitAccumulation_(const Schema& schema, bool surf)
{
  Errors::Message msg;

  if (global_op_ == Teuchos::null) {
    // constructor was given a mesh
    global_op_schema_ = schema;
    local_op_schema_ = schema;

    for (auto it = schema.begin(); it != schema.end(); ++it) {
      Teuchos::RCP<Op> op;
      // Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new
      // CompositeVectorSpace());
      // cvs->SetMesh(mesh_)->AddComponent(schema.KindToString(it->kind),
      // it->kind, it->num);
      auto cvs = CreateCompositeVectorSpace(
        mesh_, schema.KindToString(it->kind), it->kind, it->num);

      if (it->kind == AmanziMesh::CELL) {
        int old_schema = OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_CELL;
        global_op_ = Teuchos::rcp(new Operator_Cell(cvs, plist_, old_schema));
        std::string name("CELL_CELL");
        if (surf) {
          op = Teuchos::rcp(new Op_SurfaceCell_SurfaceCell(name, mesh_));
        } else {
          op = Teuchos::rcp(new Op_Cell_Cell(name, mesh_));
        }

        /*
        } else if (it->kind == AmanziMesh::FACE) {
          global_op_ = Teuchos::rcp(new Operator_Face(cvs, plist_));
          std::string name("FACE_FACE");
          op = Teuchos::rcp(new Op_Face_Face(name, mesh_));
        */

      } else if (it->kind == AmanziMesh::EDGE) {
        global_op_ = Teuchos::rcp(new Operator_Edge(cvs, plist_));
        std::string name("EDGE_EDGE");
        op = Teuchos::rcp(new Op_Edge_Edge(name, mesh_));

      } else if (it->kind == AmanziMesh::NODE) {
        global_op_ = Teuchos::rcp(new Operator_Node(cvs, plist_));
        std::string name("NODE_NODE");
        op = Teuchos::rcp(new Op_Node_Node(name, mesh_, it->num));

      } else {
        msg << "Accumulation operator: Unknown kind \""
            << schema.KindToString(it->kind) << "\".\n";
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
      int old_schema;
      Teuchos::RCP<Op> op;

      if (it->kind == AmanziMesh::CELL) {
        old_schema = OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_CELL;
        std::string name("CELL_CELL");
        if (surf) {
          op = Teuchos::rcp(new Op_SurfaceCell_SurfaceCell(name, mesh_));
        } else {
          op = Teuchos::rcp(new Op_Cell_Cell(name, mesh_));
        }

        /*
        } else if (it->kind == AmanziMesh::FACE) {
          old_schema = OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_FACE;
          std::string name("FACE_FACE");
          op = Teuchos::rcp(new Op_Face_Face(name, mesh_));
        */

      } else if (it->kind == AmanziMesh::EDGE) {
        old_schema = OPERATOR_SCHEMA_BASE_EDGE | OPERATOR_SCHEMA_DOFS_EDGE;
        std::string name("EDGE_EDGE");
        op = Teuchos::rcp(new Op_Edge_Edge(name, mesh_));

      } else if (it->kind == AmanziMesh::NODE) {
        old_schema = OPERATOR_SCHEMA_BASE_NODE | OPERATOR_SCHEMA_DOFS_NODE;
        std::string name("NODE_NODE");
        op = Teuchos::rcp(new Op_Node_Node(name, mesh_, it->num));

      } else {
        msg << "Accumulation operator: Unknown kind \""
            << schema.KindToString(it->kind) << "\".\n";
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
    mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  nfaces_owned =
    mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nnodes_owned =
    mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
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
      if (schema.base() == (*bc)->kind()) {
        Epetra_MultiVector& diag = *(*it)->diag;
        int m = diag.getNumVectors();

        for (int i = 0; i < diag.getLocalLength(); i++) {
          if (bc_model[i] == OPERATOR_BC_DIRICHLET) {
            for (int k = 0; k < m; ++k) { diag[k][i] = 0.0; }
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
    if (schema.KindToString(schema.base()) == name) return *it;
  }
  return Teuchos::null;
}

} // namespace Operators
} // namespace Amanzi
