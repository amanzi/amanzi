/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include <vector>

// TPLs
#include "Epetra_Vector.h"

// Amanzi
#include "BilinearFormFactory.hh"
#include "errors.hh"
#include "MatrixFE.hh"
#include "PreconditionerFactory.hh"
#include "WhetStoneDefs.hh"

// Amanzi::Operators
#include "Op.hh"
#include "Op_Cell_Schema.hh"
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
PDE_Elasticity::SetTensorCoefficient(
  const Teuchos::RCP<std::vector<WhetStone::Tensor>>& K)
{
  K_ = K;
  K_default_ = 1.0;
}

void
PDE_Elasticity::SetTensorCoefficient(double K)
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

  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = K_default_;

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
PDE_Elasticity::Init_(Teuchos::ParameterList& plist)
{
  // Read schema for the mimetic discretization method.
  Schema my_schema;
  Teuchos::ParameterList& schema_list = plist.sublist("schema");
  my_schema.Init(schema_list, mesh_);

  // create or check the existing Operator
  local_schema_col_ = my_schema;
  local_schema_row_ = my_schema;

  if (global_op_ == Teuchos::null) {
    global_schema_col_ = my_schema;
    global_schema_row_ = my_schema;

    // build the CVS from the global schema
    Teuchos::RCP<CompositeVectorSpace> cvs =
      Teuchos::rcp(new CompositeVectorSpace());
    cvs->SetMesh(mesh_)->SetGhosted(true);

    for (auto it = my_schema.begin(); it != my_schema.end(); ++it) {
      std::string name(my_schema.KindToString(it->kind));
      cvs->AddComponent(name, it->kind, it->num);
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

  // parse discretization  parameters
  mfd_ = WhetStone::BilinearFormFactory::Create(plist, mesh_);
}

} // namespace Operators
} // namespace Amanzi
