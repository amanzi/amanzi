/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#include <vector>

// TPLs
#include "Epetra_Vector.h"

// Amanzi
#include "BilinearFormFactory.hh"
#include "errors.hh"
#include "MatrixFE.hh"
#include "MFD3D_Elasticity.hh"
#include "WhetStoneDefs.hh"
#include "WhetStoneMeshUtils.hh"

// Amanzi::Operators
#include "Op.hh"
#include "Op_Cell_Schema.hh"
#include "Op_Node_Schema.hh"
#include "OperatorDefs.hh"
#include "Operator_Schema.hh"
#include "PDE_Elasticity.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization of the operator, scalar coefficient.
****************************************************************** */
void
PDE_Elasticity::SetTensorCoefficient(const Teuchos::RCP<std::vector<WhetStone::Tensor>>& K)
{
  K_ = K;
  K_default_.Init(1, 1);
}

void
PDE_Elasticity::SetTensorCoefficient(const WhetStone::Tensor& K)
{
  K_ = Teuchos::null;
  K_default_ = K;
}


/* ******************************************************************
* Calculate elemental matrices.
* NOTE: The input parameters are not yet used.
****************************************************************** */
void
PDE_Elasticity::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                               const Teuchos::Ptr<const CompositeVector>& p)
{
  WhetStone::DenseMatrix Acell;
  WhetStone::Tensor Kc(K_default_);

  for (int c = 0; c < ncells_owned; c++) {
    if (K_.get()) Kc = (*K_)[c];
    mfd_->StiffnessMatrix(c, Kc, Acell);
    local_op_->matrices[c] = Acell;
  }
}


/* ******************************************************************
* Put here stuff that has to be done in constructor.
****************************************************************** */
void
PDE_Elasticity::ComputeHydrostaticStress(const CompositeVector& u, CompositeVector& p)
{
  WhetStone::Tensor Tc;
  WhetStone::Tensor Kc(K_default_);
  WhetStone::DenseVector dofs;

  auto& p_c = *p.ViewComponent("cell");
  const auto& u_n = *u.ViewComponent("node", true);

  Teuchos::RCP<const Epetra_MultiVector> u_f;
  if (u.HasComponent("face")) u_f = u.ViewComponent("face", true);

  int d = mesh_->getSpaceDimension();
  auto mfd3d = Teuchos::rcp_dynamic_cast<WhetStone::MFD3D>(mfd_);

  p_c.PutScalar(0.0);
 
  for (int c = 0; c < ncells_owned; c++) {
    // nodal DoFs go first
    auto nodes = mesh_->getCellNodes(c);
    int nnodes = nodes.size();

    dofs.Reshape(d * nnodes);
    for (int n = 0; n < nnodes; ++n) {
      int v = nodes[n];
      for (int k = 0; k < d; ++k) {
        dofs(d * n + k) = u_n[k][v];
      }
    }

    // optional face DoFs
    if (u_f.get()) {
      auto faces = mesh_->getCellFaces(c);
      int nfaces = faces.size();

      dofs.Reshape(d * nnodes + nfaces);
      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        dofs(d * nnodes + n) = (*u_f)[0][f];
      }
    }

    mfd3d->H1Cell(c, dofs, Tc);

    if (K_.get()) Kc = (*K_)[c];
    WhetStone::Tensor TKc = Kc * Tc; 

    for (int k = 0; k < d; ++k) p_c[0][c] += TKc(k, k);
  }
}


/* ******************************************************************
* Cell-centered volumetric strain 
****************************************************************** */
void
PDE_Elasticity::ComputeVolumetricStrain(const CompositeVector& u, CompositeVector& e)
{
  WhetStone::Tensor Tc;
  WhetStone::DenseVector dofs;

  auto& e_c = *e.ViewComponent("cell");
  const auto& u_n = *u.ViewComponent("node", true);

  Teuchos::RCP<const Epetra_MultiVector> u_f;
  if (u.HasComponent("face")) u_f = u.ViewComponent("face", true);

  int d = mesh_->getSpaceDimension();
  auto mfd3d = Teuchos::rcp_dynamic_cast<WhetStone::MFD3D>(mfd_);

  e_c.PutScalar(0.0);
 
  for (int c = 0; c < ncells_owned; c++) {
    // nodal DoFs go first
    auto nodes = mesh_->getCellNodes(c);
    int nnodes = nodes.size();

    dofs.Reshape(d * nnodes);
    for (int n = 0; n < nnodes; ++n) {
      int v = nodes[n];
      for (int k = 0; k < d; ++k) {
        dofs(d * n + k) = u_n[k][v];
      }
    }

    // optional face DoFs
    if (u_f.get()) {
      auto faces = mesh_->getCellFaces(c);
      int nfaces = faces.size();

      dofs.Reshape(d * nnodes + nfaces);
      for (int n = 0; n < nfaces; ++n) {
        int f = faces[n];
        dofs(d * nnodes + n) = (*u_f)[0][f];
      }
    }

    mfd3d->H1Cell(c, dofs, Tc);

    for (int k = 0; k < d; ++k) e_c[0][c] += Tc(k, k);
  }
}


/* ******************************************************************
* Put here stuff that has to be done in constructor.
****************************************************************** */
void
PDE_Elasticity::Init_(Teuchos::ParameterList& plist)
{
  // generate schema for the mimetic discretization method
  Teuchos::ParameterList& schema_list = plist.sublist("schema");
  mfd_ = WhetStone::BilinearFormFactory::Create(schema_list, mesh_);

  Schema my_schema;
  base_ = my_schema.StringToKind(schema_list.get<std::string>("base"));
  my_schema.Init(mfd_, mesh_, base_);

  // create or check the existing Operator
  local_schema_col_ = my_schema;
  local_schema_row_ = my_schema;

  if (global_op_ == Teuchos::null) {
    global_schema_col_ = my_schema;
    global_schema_row_ = my_schema;

    // build the CVS from the global schema
    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true);

    for (auto it = my_schema.begin(); it != my_schema.end(); ++it) {
      int num;
      AmanziMesh::Entity_kind kind;
      std::tie(kind, std::ignore, num) = *it;

      std::string name(my_schema.KindToString(kind));
      cvs->AddComponent(name, kind, num);
    }

    global_op_ = Teuchos::rcp(new Operator_Schema(cvs, plist, my_schema));
  } else {
    // constructor was given an Operator
    global_schema_col_ = global_op_->schema_col();
    global_schema_row_ = global_op_->schema_row();
    mesh_ = global_op_->DomainMap().Mesh();
  }

  // create the local Op and register it with the global Operator
  local_op_ = Teuchos::rcp(new Op_Cell_Schema(my_schema, my_schema, mesh_));
  global_op_->OpPushBack(local_op_);

  K_ = Teuchos::null;
}


/* ******************************************************************
* Put here stuff that has to be done in constructor.
****************************************************************** */
void
PDE_Elasticity::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  // implementation of specific BCs: shear stress (2D), kinematic (3D), traction
  for (auto bc : bcs_trial_) {
    auto kind = bc->kind();
    auto type = bc->type();

    if (kind == AmanziMesh::Entity_kind::FACE &&
        type == WhetStone::DOF_Type::SCALAR) {
      ApplyBCs_ShearStress_(*bc, primary, eliminate, essential_eqn);

    } else if (kind == AmanziMesh::Entity_kind::FACE &&
               type == WhetStone::DOF_Type::POINT) {
      ApplyBCs_Traction_(*bc, primary, eliminate, essential_eqn);

    } else if (kind == AmanziMesh::NODE &&
               type == WhetStone::DOF_Type::SCALAR) {
      ApplyBCs_Kinematic_(*bc, primary, eliminate, essential_eqn);
    }
  }

  // default implementation of essential BCs
  for (auto bc : bcs_trial_) {
    auto type = bc->type();

    if (type == WhetStone::DOF_Type::POINT) {
      ApplyBCs_Cell_Point_(*bc, local_op_, primary, eliminate, essential_eqn);
    } else if (type == WhetStone::DOF_Type::SCALAR ||
               type == WhetStone::DOF_Type::NORMAL_COMPONENT) {
      ApplyBCs_Cell_Scalar_(*bc, local_op_, primary, eliminate, essential_eqn);
    }
  }
}


/* ******************************************************************
* BCs: shear stress
****************************************************************** */
void
PDE_Elasticity::ApplyBCs_ShearStress_(
  const BCs& bc, bool primary, bool eliminate, bool essential_eqn)
{
  int d = mesh_->getSpaceDimension();
  const auto& bc_model = bc.bc_model();
  const auto& bc_value = bc.bc_value();

  auto rhs = global_op_->rhs();
  Teuchos::RCP<Epetra_MultiVector> rhs_node;
  if (rhs()->HasComponent("node")) rhs_node = rhs->ViewComponent("node", true);

  rhs->PutScalarGhosted(0.0);

  for (int c = 0; c != ncells_owned; ++c) {
    auto faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];

      if (bc_model[f] == OPERATOR_BC_SHEAR_STRESS) {
        double area = mesh_->getFaceArea(f);
        const auto& tau = mesh_->getEdgeVector(f);
        auto lnodes = mesh_->getEdgeNodes(f);
        int nlnodes = lnodes.size();

        std::vector<double> weights(nlnodes, 0.5);
        if (d == 3) WhetStone::PolygonCentroidWeights(*mesh_, lnodes, area, weights);

        double value = bc_value[f];
        for (int m = 0; m < nlnodes; ++m) {
          int v = lnodes[m];
          for (int k = 0; k < d; ++k) {
            (*rhs_node)[k][v] += value * tau[k] * weights[m];
          }
        }
      }
    }
  }

  rhs->GatherGhostedToMaster(Add);
}


/* ******************************************************************
* BCs: traction
****************************************************************** */
void
PDE_Elasticity::ApplyBCs_Traction_(
  const BCs& bc, bool primary, bool eliminate, bool essential_eqn)
{
  int d = mesh_->getSpaceDimension();
  const auto& bc_model = bc.bc_model();
  const auto& bc_value = bc.bc_value_point();

  auto rhs = global_op_->rhs();
  Teuchos::RCP<Epetra_MultiVector> rhs_node;
  if (rhs()->HasComponent("node")) rhs_node = rhs->ViewComponent("node", true);

  rhs->PutScalarGhosted(0.0);

  for (int c = 0; c != ncells_owned; ++c) {
    auto faces = mesh_->getCellFaces(c);
    int nfaces = faces.size();

    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];

      if (bc_model[f] == OPERATOR_BC_NORMAL_STRESS) {
        double area = mesh_->getFaceArea(f);
        auto lnodes = mesh_->getFaceNodes(f);
        int nlnodes = lnodes.size();

        std::vector<double> weights(nlnodes, 0.5);
        if (d == 3) WhetStone::PolygonCentroidWeights(*mesh_, lnodes, area, weights);

        auto& value = bc_value[f];
        for (int m = 0; m < nlnodes; ++m) {
          int v = lnodes[m];
          for (int k = 0; k < d; ++k) {
            (*rhs_node)[k][v] += value[k] * weights[m] * area;
          }
        }
      }
    }
  }

  rhs->GatherGhostedToMaster(Add);
}


/* ******************************************************************
* BCs: kinematic
****************************************************************** */
void
PDE_Elasticity::ApplyBCs_Kinematic_(
  const BCs& bc, bool primary, bool eliminate, bool essential_eqn)
{
  int d = mesh_->getSpaceDimension();
  const auto& bc_model = bc.bc_model();
  const auto& bc_value = bc.bc_value();

  auto rhs = global_op_->rhs();
  Teuchos::RCP<Epetra_MultiVector> rhs_node;
  if (rhs()->HasComponent("node")) rhs_node = rhs->ViewComponent("node", true);

  rhs->PutScalarGhosted(0.0);

  for (int c = 0; c != ncells_owned; ++c) {
    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];
    int ncols = Acell.NumCols();

    auto nodes = mesh_->getCellNodes(c);
    int nnodes = nodes.size();

    for (int n = 0; n != nnodes; ++n) {
      int v = nodes[n];

      if (bc_model[v] == OPERATOR_BC_KINEMATIC) {
        double value = bc_value[v];
        if (local_op_->matrices_shadow[c].NumRows() == 0) {
          local_op_->matrices_shadow[c] = Acell;
        }
        auto cells = mesh_->getNodeCells(v);
        int ncells = cells.size();

        auto normal = WhetStone::getNodeUnitNormal(*mesh_, v);
        int k = (std::fabs(normal[0]) > std::fabs(normal[1])) ? 0 : 1;
        if (d == 3) k = (std::fabs(normal[k]) > std::fabs(normal[2])) ? k : 2;

        // keeps positive number on the main diagonal
        if (normal[k] < 0.0) {
          normal *= -1.0;
          value *= -1.0;
        }

        int noff(d * n);
        for (int m = 0; m < ncols; m++) Acell(noff + k, m) = 0.0;
        for (int i = 0; i < d; ++i) Acell(noff + k, noff + i) = normal[i] / ncells;
        if (v < nnodes_owned) (*rhs_node)[k][v] = value;

        if (eliminate) {
          // AMANZI_ASSERT(false);
        }
      }
    }
  }

  rhs->GatherGhostedToMaster(Add);
}

} // namespace Operators
} // namespace Amanzi
