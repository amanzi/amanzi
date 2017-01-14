/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  Base class for advection operators.
*/

#ifndef AMANZI_OPERATOR_ADVECTION_HH_
#define AMANZI_OPERATOR_ADVECTION_HH_

#include "Epetra_IntVector.h"

#include "Operator.hh"

namespace Amanzi {
namespace Operators {

class Advection {
 public:
  Advection(Teuchos::ParameterList& plist,
            Teuchos::RCP<Operator> global_op) :
      global_op_(global_op),
      mesh_(Teuchos::null) {};

  Advection(Teuchos::ParameterList& plist,
            Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      global_op_(Teuchos::null),
      mesh_(mesh) {};

  // main members
  // -- setup
  virtual void Setup(const CompositeVector& u) = 0;
  // -- populate data: arrays of local matrices
  virtual void UpdateMatrices(const CompositeVector& u) = 0;
  virtual void UpdateMatrices(const CompositeVector& u, const CompositeVector& dhdT) = 0;
  // -- boundary conditions
  virtual void ApplyBCs(const Teuchos::RCP<BCs>& bc, bool primary) = 0;

  // -- results: determine advected flux of u
  virtual void UpdateFlux(const CompositeVector& h , const CompositeVector& u,
                          const Teuchos::RCP<BCs>& bc, CompositeVector& flux) = 0;
  
  // access
  Teuchos::RCP<Operator> global_operator() { return global_op_; }
  Teuchos::RCP<const Operator> global_operator() const { return global_op_; }

 protected:
  // operators and schemas
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

