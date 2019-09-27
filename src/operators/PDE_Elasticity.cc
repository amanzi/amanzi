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
#include "BilinearFormFactory.hh"
#include "errors.hh"
#include "MatrixFE.hh"
#include "MFD3D_Elasticity.hh"
#include "PreconditionerFactory.hh"
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
void PDE_Elasticity::SetTensorCoefficient(
    const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K) {
  K_ = K;
  K_default_ = 1.0;
}

void PDE_Elasticity::SetTensorCoefficient(double K) {
  K_ = Teuchos::null;
  K_default_ = K;
}


/* ******************************************************************
* Calculate elemental matrices.
* NOTE: The input parameters are not yet used.
****************************************************************** */
void PDE_Elasticity::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                                    const Teuchos::Ptr<const CompositeVector>& p)
{
  WhetStone::DenseMatrix A, B;
  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = K_default_;
  
  // use factory for cell-based method
  if (base_ == AmanziMesh::CELL) {
    for (int c = 0; c < ncells_owned; ++c) {
    if (K_.get()) Kc = (*K_)[c];
      mfd_->StiffnessMatrix(c, Kc, A);
      local_op_->matrices[c] = A;
    }
  // special elasticity methods: there exists only one such method so far
  } else if (base_ == AmanziMesh::NODE) {
    AmanziMesh::Entity_ID_List cells;
    auto mfd_elasticity = Teuchos::rcp_dynamic_cast<WhetStone::MFD3D_Elasticity>(mfd_);

    for (int n = 0; n < nnodes_owned; ++n) {
      mesh_->node_get_cells(n, AmanziMesh::Parallel_type::ALL, &cells);

      std::vector<WhetStone::Tensor> vKc;
      for (int i = 0; i < cells.size(); ++i) {
        if (K_.get()) Kc = (*K_)[cells[i]];
        vKc.push_back(Kc);
      }
      mfd_elasticity->StiffnessMatrix_LocalStress(n, vKc, A, B);
      local_op_->matrices[n] = A;
      local_op_->matrices_shadow[n] = B;
    }
  }
}


/* ******************************************************************
* Apply boundary conditions to the local matrices. 
****************************************************************** */
void PDE_Elasticity::ApplyBCs(bool primary, bool eliminate, bool essential_eqn)
{
  if (base_ == AmanziMesh::NODE)
    ApplyBCs_Node_Point_(*bcs_trial_[0], local_op_, primary, eliminate, essential_eqn);
  else
    PDE_HelperDiscretization::ApplyBCs(primary, eliminate, essential_eqn);
}


/* ******************************************************************
* Put here stuff that has to be done in constructor.
****************************************************************** */
void PDE_Elasticity::Init_(Teuchos::ParameterList& plist)
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
  if (base_ == AmanziMesh::CELL) {
  local_op_ = Teuchos::rcp(new Op_Cell_Schema(my_schema, my_schema, mesh_));
  } else if (base_ == AmanziMesh::NODE) {
    local_op_ = Teuchos::rcp(new Op_Node_Schema(my_schema, my_schema, mesh_));
  }

  global_op_->OpPushBack(local_op_);
  
  K_ = Teuchos::null;
}


/* ******************************************************************
* Apply BCs of point type. Generalize and move to the helper class.
****************************************************************** */
void PDE_Elasticity::ApplyBCs_Node_Point_(
    const BCs& bc, Teuchos::RCP<Op> op,
    bool primary, bool eliminate, bool essential_eqn)
{
  const std::vector<int>& bc_model = bc.bc_model();
  const std::vector<std::vector<double> >& bc_value = bc.bc_value_vector();

  AmanziMesh::Entity_ID_List faces;

  CompositeVector& rhs = *global_op_->rhs();
  rhs.PutScalarGhosted(0.0);

  int d = mesh_->space_dimension(); 
  const Schema& schema_row = global_op_->schema_row();

  for (int v = 0; v != nnodes_owned; ++v) {
    WhetStone::DenseMatrix& Anode = op->matrices_shadow[v];
    int ncols = Anode.NumCols();
    int nrows = Anode.NumRows();

    mesh_->node_get_faces(v, AmanziMesh::Parallel_type::ALL, &faces);
    int nfaces = faces.size();

    // essential zero conditions for trial functions
    for (int n = 0; n != nfaces; ++n) {
      int f = faces[n];
      double area = mesh_->face_area(f) / 2;  // FIXME for 3D
      const std::vector<double>& value = bc_value[f];

      if (bc_model[f] == OPERATOR_BC_DIRICHLET) {
        for (int k = 0; k < d; ++k) {
          int noff(d*n + k);
          WhetStone::DenseVector rhs_loc(nrows);

          if (eliminate) {
            int pos = WhetStone::UniqueIndexFaceToNodes(*mesh_, f, v);
            for (int m = 0; m < nrows; m++) {
              rhs_loc(m) = Anode(m, noff) * value[d * pos + k] * area;
            }
          }

          global_op_->AssembleVectorNodeOp(v, schema_row, rhs_loc, rhs);
        }
      }
    }
  } 

  rhs.GatherGhostedToMaster(Add);
}

}  // namespace Operators
}  // namespace Amanzi

