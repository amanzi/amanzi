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

// Amanzi
#include "errors.hh"
#include "MatrixFE.hh"
#include "MFD3D_Electromagnetics.hh"
#include "Point.hh"
#include "PreconditionerFactory.hh"
#include "SuperMap.hh"
#include "WhetStoneDefs.hh"

// Amanzi::Operators
#include "PDE_Electromagnetics.hh"
#include "Op.hh"
#include "Op_Cell_Edge.hh"
#include "Op_Cell_Node.hh"
#include "OperatorDefs.hh"
#include "Operator_Edge.hh"
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

  if (local_op_schema_ == OPERATOR_SCHEMA_BASE_CELL + OPERATOR_SCHEMA_DOFS_EDGE) {
    if (K_ != Teuchos::null && K_.get()) AMANZI_ASSERT(K_->size() == ncells_owned);
  }
}


/* ******************************************************************
* Calculate elemental matrices.
* NOTE: The input parameters are not yet used.
****************************************************************** */
void PDE_Electromagnetics::UpdateMatrices()
{
  Teuchos::ParameterList plist;
  WhetStone::MFD3D_Electromagnetics mfd(plist, mesh_);
  WhetStone::DenseMatrix Acell;

  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = 1.0;
  
  for (int c = 0; c < ncells_owned; c++) {
    if (K_.get()) Kc = (*K_)[c];
    if (mfd_primary_ == WhetStone::ELECTROMAGNETICS_GENERALIZED)
      mfd.StiffnessMatrix_GradCorrection(c, Kc, Acell);
    else
      mfd.StiffnessMatrix(c, Kc, Acell);
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
  if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_CELL
                         | OPERATOR_SCHEMA_DOFS_EDGE)) {
    Teuchos::RCP<const BCs> bc_f, bc_e;
    for (auto bc = bcs_trial_.begin(); bc != bcs_trial_.end(); ++bc) {
      if ((*bc)->kind() == AmanziMesh::FACE) {
        bc_f = *bc;
      } else if ((*bc)->kind() == AmanziMesh::EDGE) {
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
  AmanziMesh::Entity_ID_List edges, faces, cells;
  std::vector<int> edirs, fdirs;

  global_op_->rhs()->PutScalarGhosted(0.0);
  Epetra_MultiVector& rhs_edge = *global_op_->rhs()->ViewComponent("edge", true);

  // support of surface integrals
  Teuchos::ParameterList plist;
  int dim = mesh_->space_dimension();
  WhetStone::MFD3D_Electromagnetics mfd3d(plist, mesh_);

  // calculate number of cells for each edge 
  // move to properties of BCs (lipnikov@lanl.gov)
  std::vector<int> edge_ncells(nedges_wghost, 0);
  for (int c = 0; c != ncells_wghost; ++c) {
    mesh_->cell_get_edges(c, &edges);
    int nedges = edges.size();

    for (int n = 0; n < nedges; ++n) {
      edge_ncells[edges[n]]++;
    }
  }

  int nn(0), nm(0);
  for (int c = 0; c != ncells_owned; ++c) {
    bool flag(true);
    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];

    // BCs of faces: typically this is magnetic flux
    if (bc_f != Teuchos::null) {
      const std::vector<int>& bc_model = bc_f->bc_model();
      const std::vector<AmanziGeometry::Point>& bc_value = bc_f->bc_value_point();

      mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);
      int nfaces = faces.size();

      for (int n = 0; n != nfaces; ++n) {
        int f = faces[n];
        const AmanziGeometry::Point& value = bc_value[f];

        if (bc_model[f] == OPERATOR_BC_NEUMANN && primary) {
          const AmanziGeometry::Point& normal = mesh_->face_normal(f);
          double area = mesh_->face_area(f);

          mesh_->face_get_edges_and_dirs(f, &edges, &edirs);
          int nedges = edges.size();

          // project magnetic flux on mesh edges
          WhetStone::DenseVector b(nedges), mb(nedges); 
          for (int i = 0; i != nedges; ++i) {
            int e = edges[i];
            const AmanziGeometry::Point& tau = mesh_->edge_vector(e);
            double len = mesh_->edge_length(e);
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

      mesh_->cell_get_edges(c, &edges);
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
  // Determine discretization
  std::string primary = plist.get<std::string>("discretization primary");
  K_symmetric_ = (plist.get<std::string>("diffusion tensor", "symmetric") == "symmetric");

  // Primary discretization methods
  if (primary == "mfd: default") {
    mfd_primary_ = WhetStone::ELECTROMAGNETICS_DEFAULT;
  } else if (primary == "mfd: generalized") {
    mfd_primary_ = WhetStone::ELECTROMAGNETICS_GENERALIZED;
  } else {
    Errors::Message msg;
    msg << "Electromagnetics: primary discretization method \"" << primary << "\" is not supported.";
    Exceptions::amanzi_throw(msg);
  }

  // Define stencil for the MFD diffusion method.
  std::vector<std::string> names;
  if (plist.isParameter("schema")) {
    names = plist.get<Teuchos::Array<std::string> > ("schema").toVector();
  } else {
    names.resize(1);
    names[0] = "edge";
    plist.set<Teuchos::Array<std::string> >("schema", names);
  }

  int dim = mesh_->space_dimension();
  int schema_dofs = 0;
  for (int i = 0; i < names.size(); i++) {
    if (names[i] == "edge" && dim == 3) {
      schema_dofs += OPERATOR_SCHEMA_DOFS_EDGE;
    } else if (names[i] == "node" && dim == 2) {
      schema_dofs += OPERATOR_SCHEMA_DOFS_NODE;
    }
  }
 
  if (schema_dofs == 0) {
    Errors::Message msg;
    msg << "Electromagnetics: \"schema\" must be EDGE (in 3D) or NODE (in 2D)";
    Exceptions::amanzi_throw(msg);
  }

  local_op_schema_ = OPERATOR_SCHEMA_BASE_CELL | schema_dofs;

  // define stencil for the assembled matrix
  int schema_prec_dofs = 0;
  if (plist.isParameter("preconditioner schema")) {
    names = plist.get<Teuchos::Array<std::string> > ("preconditioner schema").toVector();
    for (int i = 0; i < names.size(); i++) {
      if (names[i] == "edge" && dim == 3) {
        schema_prec_dofs += OPERATOR_SCHEMA_DOFS_EDGE;
      } else if (names[i] == "node" && dim == 2) {
        schema_prec_dofs += OPERATOR_SCHEMA_DOFS_NODE;
      }
    } 
  } else {
    schema_prec_dofs = schema_dofs;
  }


  // create or check the existing Operator
  int global_op_schema = schema_prec_dofs;  
  if (global_op_ == Teuchos::null) {
    global_op_schema_ = global_op_schema;

    // build the CVS from the global schema
    Teuchos::RCP<CompositeVectorSpace> cvs = Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true);

    if (global_op_schema & OPERATOR_SCHEMA_DOFS_EDGE) {
      cvs->AddComponent("edge", AmanziMesh::EDGE, 1);
    } else if (global_op_schema & OPERATOR_SCHEMA_DOFS_NODE) {
      cvs->AddComponent("node", AmanziMesh::NODE, 1);
    }

    // choose the Operator from the prec schema
    Teuchos::ParameterList operator_list = plist.sublist("operator");
    if (schema_prec_dofs == OPERATOR_SCHEMA_DOFS_EDGE) {
      global_op_ = Teuchos::rcp(new Operator_Edge(cvs, plist));
    } else if (schema_prec_dofs == OPERATOR_SCHEMA_DOFS_NODE) {
      global_op_ = Teuchos::rcp(new Operator_Node(cvs, plist));
    } else {
      Errors::Message msg;
      msg << "Electromagnetics: \"preconditioner schema\" must be EDGE";
      Exceptions::amanzi_throw(msg);
    }

  } else {
    // constructor was given an Operator
    global_op_schema_ = global_op_->schema();
    mesh_ = global_op_->DomainMap().Mesh();
  }

  // create the local Op and register it with the global Operator
  if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_EDGE)) {
    std::string name = "Electromagnetics: CELL_EDGE";
    local_op_ = Teuchos::rcp(new Op_Cell_Edge(name, mesh_));
  } else if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_NODE)) {
    std::string name = "Electromagnetics: CELL_NODE";
    local_op_ = Teuchos::rcp(new Op_Cell_Node(name, mesh_));
  } else {
    AMANZI_ASSERT(0);
  }
  global_op_->OpPushBack(local_op_);

  K_ = Teuchos::null;
}

}  // namespace Operators
}  // namespace Amanzi
