/*
  This is the operators component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  Discrete advection operator of a surface.
*/

#ifndef AMANZI_OPERATOR_ADVECTION_HH_
#define AMANZI_OPERATOR_ADVECTION_HH_

#include "Epetra_IntVector.h"

#include "Operator.hh"


namespace Amanzi {
namespace Operators {

class OperatorAdvection {
 public:

  OperatorAdvection(Teuchos::ParameterList& plist,
                    Teuchos::RCP<Operator> global_op) :
      global_op_(global_op),
      mesh_(Teuchos::null)
  {
    InitAdvection_(plist);
  }

  OperatorAdvection(Teuchos::ParameterList& plist,
                    Teuchos::RCP<AmanziMesh::Mesh> mesh) :
      global_op_(Teuchos::null),
      mesh_(mesh)
  {
    InitAdvection_(plist);
  }

  OperatorAdvection(Teuchos::ParameterList& plist,
                    Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      global_op_(Teuchos::null),
      mesh_(mesh)
  {
    InitAdvection_(plist);
  }

  // main members
  void Setup(const CompositeVector& u);
  void UpdateMatrices(const CompositeVector& u);
  void UpdateMatrices(const CompositeVector& u,
                      const CompositeVector& dhdT);

 protected:
  void InitAdvection_(Teuchos::ParameterList& plist);
  void IdentifyUpwindCells_(const CompositeVector& u);

 protected:
  Teuchos::RCP<Epetra_IntVector> upwind_cell_, downwind_cell_;

  // operator
  Teuchos::RCP<Operator> global_op_;
  Teuchos::RCP<Op> local_op_;
  int global_op_schema_, local_op_schema_;

  // mesh info
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;
  int nnodes_owned, nnodes_wghost;
  
};

}  // namespace Operators
}  // namespace Amanzi

#endif

