/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  This operator is a collection of local "DIAGONAL" Ops.
*/

//#include "WhetStoneMeshUtils.hh"

#include "OperatorUtils.hh"
#include "Operator_Cell.hh"
//#include "Operator_Edge.hh"
//#include "Operator_Node.hh"
#include "Op_Cell_Cell.hh"
//#include "Op_Edge_Edge.hh"
//#include "Op_Node_Node.hh"
#include "Op_SurfaceCell_SurfaceCell.hh"
#include "PDE_Accumulation.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Modifier for diagonal operators.  Op += du
****************************************************************** */
void PDE_Accumulation::AddAccumulationTerm(
    const CompositeVector& du, const std::string& name)
{
  Teuchos::RCP<Op> op = FindOp_(name);
  op->diag->update(1.,*du.GetComponent(name, false),1.);
}


/* ******************************************************************
* Modifier for diagonal operators.  Op += du * vol / dt.
****************************************************************** */
void PDE_Accumulation::AddAccumulationTerm(
    const CompositeVector& du, double dT, const std::string& name, bool volume)
{
  Teuchos::RCP<Op> op = FindOp_(name);
  auto& diag = op->diag;
  auto du_c = du.GetComponent(name, false);

  int n = du_c->getLocalLength();
  int m = du_c->getNumVectors();

  if (volume) {
    MultiVector_type vol(du_c->getMap(), 1);
    assert(false); 
    //CalculateEntityVolume_(vol, name);
    //op->diag->elementWiseMultiply(1.0/dT, vol, *du.GetComponent(name, false), 1.0);

  } else {
    op->diag->update(1.0/dT, *du.GetComponent(name, false), 1.0);
  }
}


/* ******************************************************************
* Modifier for diagonal operators and rhs.
* Op  += alpha * s1 * vol
* Rhs += alpha * s2 * vol
****************************************************************** */
void PDE_Accumulation::AddAccumulationRhs(
    const CompositeVector& s1,
    const CompositeVector& s2,
    double alpha,
    const std::string& name,
    bool volume)
{
  Teuchos::RCP<Op> op = FindOp_(name);
  MultiVector_type& diag = *op->diag;

  const auto& s1c = s1.ViewComponent(name);
  const auto& s2c = s2.ViewComponent(name);  

  int n = s1c.size();
  int m = s1c.extent(1);

  auto rhs = global_operator()->rhs()->ViewComponent(name);

  AMANZI_ASSERT(s1c.size() == s2c.size());
  AMANZI_ASSERT(s1c.size() == diag.getLocalLength());  
  AMANZI_ASSERT(s2c.size() == rhs.size());

  AMANZI_ASSERT(s1c.size() == s2c.size());
  AMANZI_ASSERT(s1c.extent(1) == diag.getNumVectors());  
  AMANZI_ASSERT(s2c.size() == rhs.size());
  

  if (volume) {
    CompositeVector vol(s1);
    CalculateEntityVolume_(vol, name);
    auto volc = vol.ViewComponent(name); 

    for (int k = 0; k < m; k++) {
      auto diag_data = diag.getDataNonConst(k); 
      for (int i = 0; i < n; i++) {
        diag_data[i] += volc(i,0) * s1c(k,i) * alpha;
        rhs(k,i) += volc(i,0) * s2c(k,i) * alpha;
      } 
    }

  } else {
    for (int k = 0; k < m; k++) {
      auto diag_data = diag.getDataNonConst(k); 
      for (int i = 0; i < n; i++) {
        diag_data[i] += s1c(k,i) * alpha;
        rhs(k,i) += s2c(k,i) * alpha;        
      } 
    }
  }
}


/* ******************************************************************
* Linearized update methods with storage terms for component "name".
* Op  += ss * vol / dt
* RHS += s0 * vol * u0 / dt
****************************************************************** */
void PDE_Accumulation::AddAccumulationDelta(
    const CompositeVector& u0,
    const CompositeVector& s0, const CompositeVector& ss,
    double dT, const std::string& name)
{
  Teuchos::RCP<Op> op = FindOp_(name);
  MultiVector_type& diag = *op->diag;

  CompositeVector vol(ss);
  CalculateEntityVolume_(vol, name);

  const auto& u0c = u0.ViewComponent(name);
  const auto& s0c = s0.ViewComponent(name);
  const auto& ssc = ss.ViewComponent(name);

  auto volc = vol.ViewComponent(name); 
  auto rhs = global_operator()->rhs()->ViewComponent(name);
  auto diag_view = diag.getLocalViewHost(); 


  int n = u0c.size();
  int m = u0c.extent(1); 
  for (int k = 0; k < m; ++k) {
    for (int i = 0; i < n; i++) {
      double factor = volc(i,0) / dT;
      diag_view(i,k) += factor * ssc(i,k);
      rhs(i,k) += factor * s0c(i,k) * u0c(i,k);
    }
  }
}


/* ******************************************************************
* Linearized update methods with storage terms for component "name".
* Op  += vol / dt
* RHS += vol * u0 / dt
****************************************************************** */
void PDE_Accumulation::AddAccumulationDelta(
    const CompositeVector& u0,
    double dT, const std::string& name)
{
  Teuchos::RCP<Op> op = FindOp_(name);
  MultiVector_type& diag = *op->diag;

  CompositeVector vol(u0);
  CalculateEntityVolume_(vol, name);

  const auto& u0c = u0.ViewComponent(name);
  auto volc = vol.ViewComponent(name); 
  auto rhs = global_operator()->rhs()->ViewComponent(name);

  int n = u0c.size();
  int m = u0c.extent(1);
  for (int k = 0; k < m; ++k) {
    auto diag_data = diag.getDataNonConst(k); 
    for (int i = 0; i < n; i++) {
      double factor = volc(0,i) / dT;
      diag_data[i] += factor;
      rhs(k,i) += factor * u0c(i,k);
    }
  }
}


/* ******************************************************************
* Linearized update methods with storage terms for component "name".
* Op  += ss
* RHS += ss * u0
****************************************************************** */
void PDE_Accumulation::AddAccumulationDeltaNoVolume(
    const CompositeVector& u0, const CompositeVector& ss, const std::string& name)
{
  if (!ss.HasComponent(name)) AMANZI_ASSERT(false);

  Teuchos::RCP<Op> op = FindOp_(name);
  MultiVector_type& diag = *op->diag;

  const auto& u0c = u0.ViewComponent(name);
  const auto& ssc = ss.ViewComponent(name);

  auto rhs = global_operator()->rhs()->ViewComponent(name);

  int n = u0c.size();
  int m = u0c.extent(1);
  for (int k = 0; k < m; ++k) {
    auto diag_data = diag.getDataNonConst(k); 
    for (int i = 0; i < n; i++) {
      diag_data[i] += ssc(k,i);
      rhs(k,i) += ssc(k,i) * u0c(k,i);
    }
  }
}


/* ******************************************************************
* Calculate entity volume for component "name" of vector ss.
****************************************************************** */
void PDE_Accumulation::CalculateEntityVolume_(
    CompositeVector& volume, const std::string& name)
{
  AmanziMesh::Entity_ID_List nodes;

  if (name == "cell" && volume.HasComponent("cell")) {
    auto vol = volume.ViewComponent(name); 

    for (int c = 0; c != ncells_owned; ++c) {
      vol(c,0) = mesh_->cell_volume(c); 
    }

  } else if (name == "face" && volume.HasComponent("face")) {
    // Missing code.
    AMANZI_ASSERT(false);

  } else if (name == "edge" && volume.HasComponent("edge")) {
    auto vol = volume.ViewComponent(name, true); 
    //vol.putScalar(0.0);
    for(int c = 0 ; c != ncells_owned; ++c){
      vol(0,c) = 0.0; 
    }

    for (int c = 0; c != ncells_owned; ++c) {
      Kokkos::View<AmanziMesh::Entity_ID*> edges;  
      mesh_->cell_get_edges(c, edges);
      int nedges = edges.size();

      for (int i = 0; i < nedges; i++) {
        vol(0,edges[i]) += mesh_->cell_volume(c) / nedges; 
      }
    }
    volume.GatherGhostedToMaster(name);

  } else if (name == "node" && volume.HasComponent("node")) {
    AMANZI_ASSERT(0 && "PDE_Accumulation on nodes not yet implemented");
    // MultiVector_type& vol = *volume.ViewComponent(name, true); 
    // vol.PutScalar(0.0);

    // for (int c = 0; c != ncells_owned; ++c) {
    //   mesh_->cell_get_nodes(c, &nodes);
    //   int nnodes = nodes.size();

    //   double cellvolume = mesh_->cell_volume(c);
    //   std::vector<double> weights(nnodes, 1.0 / nnodes);

    //   if (mesh_->space_dimension() == 2) {
    //     WhetStone::PolygonCentroidWeights(*mesh_, nodes, cellvolume, weights);
    //   }

    //   for (int i = 0; i < nnodes; i++) {
    //     vol[0][nodes[i]] += weights[i] * cellvolume; 
    //   }
    // }
    // volume.GatherGhostedToMaster(name);

  } else {
    AMANZI_ASSERT(false);
  }
}


/* ******************************************************************
* Note: When complex schema is used to create a set of local ops, the
* the local local_op_ is not well defined.
****************************************************************** */
void PDE_Accumulation::InitAccumulation_(const Schema& schema, bool surf)
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
      Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
      cvs->SetMesh(mesh_)->SetGhosted(true);
      cvs->AddComponent(schema.KindToString(kind), kind, num);
      
      if (kind == AmanziMesh::CELL) {
        int old_schema = OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_CELL;
        global_op_ = Teuchos::rcp(new Operator_Cell(cvs->CreateSpace(), plist_, old_schema));
        std::string name("CELL_CELL");
        if (surf) {
          op = Teuchos::rcp(new Op_SurfaceCell_SurfaceCell(name, mesh_));
        } else {
          op = Teuchos::rcp(new Op_Cell_Cell(name, mesh_));
        }

      /*
      } else if (kind == AmanziMesh::FACE) {
        global_op_ = Teuchos::rcp(new Operator_Face(cvs, plist_));
        std::string name("FACE_FACE");
        op = Teuchos::rcp(new Op_Face_Face(name, mesh_));
      */

      } else if (kind == AmanziMesh::EDGE) {
        assert(false); 
        //global_op_ = Teuchos::rcp(new Operator_Edge(cvs, plist_));
        //std::string name("EDGE_EDGE");
        //op = Teuchos::rcp(new Op_Edge_Edge(name, mesh_));

      } else if (kind == AmanziMesh::NODE) {
        assert(false); 
        //global_op_ = Teuchos::rcp(new Operator_Node(cvs, plist_));
        //std::string name("NODE_NODE");
        //op = Teuchos::rcp(new Op_Node_Node(name, mesh_, num));

      } else {
        msg << "Accumulation operator: Unknown kind \"" << schema.KindToString(kind) << "\".\n";
        Exceptions::amanzi_throw(msg);
      }

      global_op_->OpPushBack(op);
      local_ops_.push_back(op);
    }

  } else {
    // constructor was given an Operator
    global_op_schema_ = global_op_->schema_row();
    mesh_ = global_op_->getDomainMap()->Mesh();

    for (auto it = schema.begin(); it != schema.end(); ++it) {
      std::tie(kind, std::ignore, num) = *it;

      int old_schema;
      Teuchos::RCP<Op> op;

      if (kind == AmanziMesh::CELL) {
        old_schema = OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_CELL;
        std::string name("CELL_CELL");
        if (surf) {
          op = Teuchos::rcp(new Op_SurfaceCell_SurfaceCell(name, mesh_));
        } else {
          op = Teuchos::rcp(new Op_Cell_Cell(name, mesh_));
        }

      /*
      } else if (kind == AmanziMesh::FACE) {
        old_schema = OPERATOR_SCHEMA_BASE_FACE | OPERATOR_SCHEMA_DOFS_FACE;
        std::string name("FACE_FACE");
        op = Teuchos::rcp(new Op_Face_Face(name, mesh_));
      */

      } else if (kind == AmanziMesh::EDGE) {
        assert(false); 
        //old_schema = OPERATOR_SCHEMA_BASE_EDGE | OPERATOR_SCHEMA_DOFS_EDGE;
        //std::string name("EDGE_EDGE");
        //op = Teuchos::rcp(new Op_Edge_Edge(name, mesh_));

      } else if (kind == AmanziMesh::NODE) {
        assert(false); 
        //old_schema = OPERATOR_SCHEMA_BASE_NODE | OPERATOR_SCHEMA_DOFS_NODE;
        //std::string name("NODE_NODE");
        //op = Teuchos::rcp(new Op_Node_Node(name, mesh_, num));

      } else {
        msg << "Accumulation operator: Unknown kind \"" << schema.KindToString(kind) << "\".\n";
        Exceptions::amanzi_throw(msg);
      }

      // register the accumulation Op
      local_op_schema_.Init(old_schema);
      global_op_->OpPushBack(op);
      local_ops_.push_back(op);
    }
  }

  // mesh info
  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
}


/* ******************************************************************
* Apply boundary conditions to 
****************************************************************** */
void PDE_Accumulation::ApplyBCs()
{
  for (auto bc = bcs_trial_.begin(); bc != bcs_trial_.end(); ++bc) {
    const auto bc_model = (*bc)->bc_model();

    for (auto it = local_ops_.begin(); it != local_ops_.end(); ++it) {
      const Schema& schema = local_op_schema_;//(*it)->schema_row();
      if (schema.base() == (*bc)->kind()) {
        MultiVector_type& diag = *(*it)->diag;
        int m = diag.getNumVectors();
        for (int i = 0; i < diag.getLocalLength(); i++) {
          if (bc_model[i] == OPERATOR_BC_DIRICHLET) {
            auto diag_data = diag.getDataNonConst(i); 
            for (int k = 0; k < m; ++k) {
              diag_data[k] = 0.0;
            }
          }
        }
      }
    }
  }
}


/* ******************************************************************
* Return operator with the given base.
****************************************************************** */
Teuchos::RCP<Op> PDE_Accumulation::FindOp_(const std::string& name) const
{
  for (auto it = local_ops_.begin(); it != local_ops_.end(); ++it) {
    const Schema& schema = local_op_schema_;//(*it)->schema_row();
    if (schema.KindToString(schema.base()) == name) 
      return *it;
  }
  return Teuchos::null;
}

}  // namespace Operators
}  // namespace Amanzi



