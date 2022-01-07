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
#include "MFD3D_Electromagnetics.hh"
#include "Point.hh"
#include "SuperMap.hh"
#include "WhetStoneDefs.hh"

// Amanzi::Operators
#include "PDE_Electromagnetics.hh"
#include "Op.hh"
#include "Op_Cell_Edge.hh"
#include "Op_Cell_Node.hh"
#include "OperatorDefs.hh"
#include "Operator_Edge.hh"
#include "Operator_Factory.hh"
#include "Operator_Node.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization of the operator, scalar coefficient.
****************************************************************** */
void PDE_Electromagnetics::SetTensorCoefficient(
    const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K)
{
  K_ = K;
  if (K_ != Teuchos::null && K_.get() && mesh_->getSpaceDimension() == 3) 
    AMANZI_ASSERT(K_->size() == ncells_owned);
}


/* ******************************************************************
* Calculate elemental matrices.
* NOTE: The input parameters are not yet used.
****************************************************************** */
void PDE_Electromagnetics::UpdateMatrices(
    const Teuchos::Ptr<const CompositeVector>& u,
    const Teuchos::Ptr<const CompositeVector>& p)
{
  WhetStone::DenseMatrix Acell;

  WhetStone::Tensor Kc(mesh_->getSpaceDimension(), 1);
  Kc(0, 0) = 1.0;
  
  for (int c = 0; c < ncells_owned; c++) {
    if (K_.get()) Kc = (*K_)[c];
    mfd_->StiffnessMatrix(c, Kc, Acell);
    local_op_->matrices[c] = Acell;
  }
}


/* ******************************************************************
* Apply boundary conditions to the local matrices. We always zero-out
* matrix rows for essential test BCs. As to trial BCs, there are
* options: (a) eliminate or not, (b) if eliminate, then put 1 on
* the diagonal or not.
****************************************************************** */
void PDE_Electromagnetics::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  if (local_schema_col_.get_base() == AmanziMesh::Entity_kind::CELL && mesh_->getSpaceDimension() == 3) {
    Teuchos::RCP<const BCs> bc_f, bc_e;
    for (auto bc = bcs_trial_.begin(); bc != bcs_trial_.end(); ++bc) {
      if ((*bc)->kind() == AmanziMesh::Entity_kind::FACE) {
        bc_f = *bc;
      } else if ((*bc)->kind() == AmanziMesh::Entity_kind::EDGE) {
        bc_e = *bc;
      }
    }
    ApplyBCs_Edge_(bc_f.ptr(), bc_e.ptr(), primary, eliminate, essential_eqn);
  }
}


/* ******************************************************************
* Apply BCs on cell operators
****************************************************************** */
void PDE_Electromagnetics::ApplyBCs_Edge_(
    const Teuchos::Ptr<const BCs>& bc_f,
    const Teuchos::Ptr<const BCs>& bc_e,
    bool primary, bool eliminate, bool essential_eqn)
{
  AmanziMesh::Entity_ID_List edges, cells;
  std::vector<int> edirs;

  global_op_->rhs()->PutScalarGhosted(0.0);
  Epetra_MultiVector& rhs_edge = *global_op_->rhs()->ViewComponent("edge", true);

  // support of surface integrals
  Teuchos::ParameterList plist;
  int dim = mesh_->getSpaceDimension();
  WhetStone::MFD3D_Electromagnetics mfd3d(plist, mesh_);

  // calculate number of cells for each edge 
  // move to properties of BCs (lipnikov@lanl.gov)
  std::vector<int> edge_ncells(nedges_wghost, 0);
  for (int c = 0; c != ncells_wghost; ++c) {
    mesh_->getCellEdges(c, edges);
    int nedges = edges.size();

    for (int n = 0; n < nedges; ++n) {
      edge_ncells[edges[n]]++;
    }
  }

  for (int c = 0; c != ncells_owned; ++c) {
    bool flag(true);
    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];

    // BCs of faces: typically this is magnetic flux
    if (bc_f != Teuchos::null) {
      const std::vector<int>& bc_model = bc_f->bc_model();
      const std::vector<AmanziGeometry::Point>& bc_value = bc_f->bc_value_point();

      const auto& faces = mesh_->getCellFaces(c);
      const auto& fdirs = mesh_->getCellFaceDirections(c);
      int nfaces = faces.size();

      for (int n = 0; n != nfaces; ++n) {
        int f = faces[n];
        const AmanziGeometry::Point& value = bc_value[f];

        if (bc_model[f] == OPERATOR_BC_NEUMANN && primary) {
          const AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
          double area = mesh_->getFaceArea(f);

          mesh_->getFaceEdgesAndDirs(f, edges, &edirs);
          int nedges = edges.size();

          // project magnetic flux on mesh edges
          WhetStone::DenseVector b(nedges), mb(nedges); 
          for (int i = 0; i != nedges; ++i) {
            int e = edges[i];
            const AmanziGeometry::Point& tau = mesh_->getEdgeVector(e);
            double len = mesh_->getEdgeLength(e);
            b(i) = ((value^normal) * tau) / (area * len) * edirs[i];
          }

          // calculate inner product matrix
          WhetStone::Tensor T(dim, 1);
          T(0, 0) = 1.0;

          WhetStone::DenseMatrix M(nedges, nedges);
          mfd3d.MassMatrixBoundary(f, T, M);
          M.Multiply(b, mb, false);

          // assemble data in the right-hand side
          for (int i = 0; i != nedges; ++i) {
            int e = edges[i];
            rhs_edge[0][e] -= mb(i) * edirs[i] * fdirs[n];
          }
        }
      }
    }

    // BCs of edges: typically this is electric field
    if (bc_e != Teuchos::null) {
      const std::vector<int>& bc_model = bc_e->bc_model();
      const std::vector<double>& bc_value = bc_e->bc_value();

      mesh_->getCellEdges(c, edges);
      int nedges = edges.size();

      // essential conditions for test functions
      for (int n = 0; n != nedges; ++n) {
        int e = edges[n];
        if (bc_model[e] == OPERATOR_BC_DIRICHLET) {
          if (flag) {  // make a copy of elemental matrix
            local_op_->matrices_shadow[c] = Acell;
            flag = false;
          }
          for (int m = 0; m < nedges; m++) Acell(n, m) = 0.0;
        }
      }

      for (int n = 0; n != nedges; ++n) {
        int e = edges[n];
        double value = bc_value[e];

        if (bc_model[e] == OPERATOR_BC_DIRICHLET) {
          if (flag) {  // make a copy of cell-based matrix
            local_op_->matrices_shadow[c] = Acell;
            flag = false;
          }
     
          if (eliminate) {
            for (int m = 0; m < nedges; m++) {
              rhs_edge[0][edges[m]] -= Acell(m, n) * value;
              Acell(m, n) = 0.0;
            }
          }

          if (essential_eqn) {
            if (e < nedges_owned) rhs_edge[0][e] = value;
            Acell(n, n) = 1.0 / edge_ncells[e];
          }
        }
      }
    }
  } 

  global_op_->rhs()->GatherGhostedToMaster("edge", Add);
}


/* ******************************************************************
* Put here stuff that has to be done in constructor.
****************************************************************** */
void PDE_Electromagnetics::Init_(Teuchos::ParameterList& plist)
{
  int dim = mesh_->getSpaceDimension();

  // domain and range of this operator are equal
  Teuchos::ParameterList domain = plist.sublist("schema electric");
  auto base = global_schema_row_.StringToKind(domain.get<std::string>("base"));

  // discretization method:
  mfd_ = WhetStone::BilinearFormFactory::Create(domain, mesh_);

  if (global_op_ == Teuchos::null) {
    // constructor was given a mesh
    local_schema_col_.Init(mfd_, mesh_, base);
    global_schema_col_ = local_schema_col_;

    local_schema_row_ = local_schema_col_;
    global_schema_row_ = global_schema_col_;

    Operator_Factory factory;
    factory.set_mesh(mesh_);
    factory.set_plist(Teuchos::rcpFromRef(plist));
    factory.set_schema(global_schema_row_);

    global_op_ = factory.CreateFromSchema();

  } else {
    // constructor was given an Operator
    global_schema_row_ = global_op_->schema_row();
    global_schema_col_ = global_op_->schema_col();

    mesh_ = global_op_->DomainMap().Mesh();
    local_schema_col_.Init(mfd_, mesh_, base);
    local_schema_row_ = local_schema_col_;
  }

  // create the local Op and register it with the global Operator
  if (local_schema_col_.get_base() == AmanziMesh::Entity_kind::CELL && dim == 3) {
    std::string name = "Electromagnetics: CELL_EDGE";
    local_op_ = Teuchos::rcp(new Op_Cell_Edge(name, mesh_));
  } else if (local_schema_col_.get_base() == AmanziMesh::Entity_kind::CELL && dim == 2) {
    std::string name = "Electromagnetics: CELL_NODE";
    local_op_ = Teuchos::rcp(new Op_Cell_Node(name, mesh_));
  } else {
    AMANZI_ASSERT(0);
  }
  global_op_->OpPushBack(local_op_);

  // other parameters
  K_ = Teuchos::null;
  K_symmetric_ = (plist.get<std::string>("diffusion tensor", "symmetric") == "symmetric");
}


/* ******************************************************************
* Additional data for AMS solver: coordinates of nodes
****************************************************************** */
Teuchos::RCP<Epetra_MultiVector> PDE_Electromagnetics::GraphGeometry()
{
  int d = mesh_->getSpaceDimension();
  auto map = mesh_->getMap(AmanziMesh::Entity_kind::NODE, false);
  auto xyz = Teuchos::rcp(new Epetra_MultiVector(map, d));

  AmanziGeometry::Point xv;
  for (int n = 0; n < nnodes_owned; ++n) {
    xv = mesh_->getNodeCoordinate(n);
    for (int i = 0; i < d; ++i) (*xyz)[i][n] = xv[i];
  }

  return xyz;
}


/* ******************************************************************
* Additional data for AMS solver: gradient
****************************************************************** */
Teuchos::RCP<Epetra_CrsMatrix> PDE_Electromagnetics::GradientOperator()
{
  auto map_row = mesh_->getMap(AmanziMesh::Entity_kind::EDGE, false);
  auto map_col = mesh_->getMap(AmanziMesh::Entity_kind::NODE, false);
  auto map_col_wghost = mesh_->getMap(AmanziMesh::Entity_kind::NODE, true);
  auto G = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map_row, map_col_wghost, 2));

  int ierr(0), n1, n2, lid_c[2];
  double values[2];

  for (int e = 0; e != nedges_owned; ++e) {
    double len = mesh_->getEdgeLength(e);
    mesh_->getEdgeNodes(e, &n1, &n2);

    lid_c[0] = map_col_wghost.GID(n1); 
    lid_c[1] = map_col_wghost.GID(n2); 

    values[0] =-1.0 / len;
    values[1] = 1.0 / len;

    ierr |= G->InsertGlobalValues(map_row.GID(e), 2, values, lid_c);
  }
  AMANZI_ASSERT(!ierr);

  G->FillComplete(map_col, map_row);
  return G;
}

}  // namespace Operators
}  // namespace Amanzi
