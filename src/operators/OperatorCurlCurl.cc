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
#include "LinearOperator.hh"
#include "LinearOperatorFactory.hh"
#include "MatrixFE.hh"
#include "mfd3d_electromagnetics.hh"
#include "PreconditionerFactory.hh"
#include "SuperMap.hh"
#include "WhetStoneDefs.hh"

// Operators
#include "Op.hh"
#include "Op_Cell_Edge.hh"
#include "OperatorDefs.hh"
#include "Operator_Edge.hh"

#include "OperatorCurlCurl.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Initialization of the operator, scalar coefficient.
****************************************************************** */
void OperatorCurlCurl::SetTensorCoefficient(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K)
{
  K_ = K;

  if (local_op_schema_ == OPERATOR_SCHEMA_BASE_CELL + OPERATOR_SCHEMA_DOFS_EDGE) {
    if (K_ != Teuchos::null && K_.get()) ASSERT(K_->size() == ncells_owned);
  }
}


/* ******************************************************************
* Calculate elemental matrices.
****************************************************************** */
void OperatorCurlCurl::UpdateMatrices()
{
  AmanziMesh::Entity_ID_List edges;
  WhetStone::MFD3D_Electromagnetics mfd(mesh_);

  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = 1.0;
  
  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_edges(c, &edges);
    int nedges = edges.size();

    WhetStone::DenseMatrix Acell(nedges, nedges);
    if (K_.get()) Kc = (*K_)[c];
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
void OperatorCurlCurl::ApplyBCs(bool primary, bool eliminate)
{
  if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_CELL
                         | OPERATOR_SCHEMA_DOFS_EDGE)) {
    Teuchos::RCP<BCs> bc_f, bc_e;
    for (std::vector<Teuchos::RCP<BCs> >::iterator bc = bcs_trial_.begin();
        bc != bcs_trial_.end(); ++bc) {
      if ((*bc)->type() == OPERATOR_BC_TYPE_FACE) {
        bc_f = *bc;
      } else if ((*bc)->type() == OPERATOR_BC_TYPE_EDGE) {
        bc_e = *bc;
      }
    }
    ApplyBCs_Edge_(bc_f.ptr(), bc_e.ptr(), primary, eliminate);
  }
}


/* ******************************************************************
* Apply BCs on cell operators
****************************************************************** */
void OperatorCurlCurl::ApplyBCs_Edge_(const Teuchos::Ptr<BCs>& bc_f,
                                      const Teuchos::Ptr<BCs>& bc_e,
                                      bool primary, bool eliminate)
{
  AmanziMesh::Entity_ID_List edges, cells;

  global_op_->rhs()->PutScalarGhosted(0.0);
  Epetra_MultiVector& rhs_edge = *global_op_->rhs()->ViewComponent("edge", true);

  int nn(0), nm(0);
  for (int c = 0; c != ncells_owned; ++c) {
    bool flag(true);
    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];

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

          if (primary) {
            // mesh_->edge_get_cells(e, AmanziMesh::USED, &cells);
            int ncells = 4;
            rhs_edge[0][e] += value / ncells;
            Acell(n, n) = 1.0 / ncells;
          }
        }
      }
    }
  } 

  global_op_->rhs()->GatherGhostedToMaster("edge", Add);
}


/* ******************************************************************
* TBW 
* **************************************************************** */
void OperatorCurlCurl::UpdateFields(const CompositeVector& E, CompositeVector& B)
{
  B.PutScalar(0.0);
}


/* ******************************************************************
* Put here stuff that has to be done in constructor.
****************************************************************** */
void OperatorCurlCurl::InitCurlCurl_(Teuchos::ParameterList& plist)
{
  // Determine discretization
  std::string primary = plist.get<std::string>("discretization primary");
  K_symmetric_ = (plist.get<std::string>("diffusion tensor", "symmetric") == "symmetric");

  // Primary discretization methods
  if (primary == "mfd: optimized for sparsity") {
    mfd_primary_ = WhetStone::DIFFUSION_OPTIMIZED_FOR_SPARSITY;
  } else if (primary == "mfd: default") {
    mfd_primary_ = WhetStone::DIFFUSION_POLYHEDRA_SCALED;
  } else {
    Errors::Message msg;
    msg << "OperatorCurlCurl: primary discretization method \"" << primary << "\" is not supported.";
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

  int schema_dofs = 0;
  for (int i = 0; i < names.size(); i++) {
    if (names[i] == "edge") {
      schema_dofs += OPERATOR_SCHEMA_DOFS_EDGE;
    }
  }

  if (schema_dofs == OPERATOR_SCHEMA_DOFS_EDGE) {
    local_op_schema_ = OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_EDGE;
  } else {
    Errors::Message msg;
    msg << "OperatorDiffusion: \"schema\" must be CELL, FACE+CELL, or NODE";
    Exceptions::amanzi_throw(msg);
  }

  // define stencil for the assembled matrix
  int schema_prec_dofs = 0;
  if (plist.isParameter("preconditioner schema")) {
    names = plist.get<Teuchos::Array<std::string> > ("preconditioner schema").toVector();
    for (int i = 0; i < names.size(); i++) {
      if (names[i] == "edge") {
        schema_prec_dofs += OPERATOR_SCHEMA_DOFS_EDGE;
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

    if (global_op_schema & OPERATOR_SCHEMA_DOFS_EDGE)
      cvs->AddComponent("edge", AmanziMesh::EDGE, 1);

    // choose the Operator from the prec schema
    Teuchos::ParameterList operator_list = plist.sublist("operator");
    if (schema_prec_dofs == OPERATOR_SCHEMA_DOFS_EDGE) {
      global_op_ = Teuchos::rcp(new Operator_Edge(cvs, plist));
    } else {
      Errors::Message msg;
      msg << "OperatorCurlCurl: \"preconditioner schema\" must be EDGE";
      Exceptions::amanzi_throw(msg);
    }

  } else {
    // constructor was given an Operator
    global_op_schema_ = global_op_->schema();
    mesh_ = global_op_->DomainMap().Mesh();
  }

  // create the local Op and register it with the global Operator
  if (local_op_schema_ == (OPERATOR_SCHEMA_BASE_CELL | OPERATOR_SCHEMA_DOFS_EDGE)) {
      std::string name = "CurlCurl: CELL_EDGE";
      local_op_ = Teuchos::rcp(new Op_Cell_Edge(name, mesh_));
  } else {
    ASSERT(0);
  }
  global_op_->OpPushBack(local_op_);
  
  // mesh info
  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nedges_owned = mesh_->num_entities(AmanziMesh::EDGE, AmanziMesh::OWNED);

  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  nedges_wghost = mesh_->num_entities(AmanziMesh::EDGE, AmanziMesh::USED);

  K_ = Teuchos::null;
}

}  // namespace Operators
}  // namespace Amanzi
