/*
  This is the Operator component of the Amanzi code.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
          Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_ACCUMULATION_HH_
#define AMANZI_OPERATOR_ACCUMULATION_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "CompositeVector.hh"

#include "Operator.hh"
#include "OperatorDefs.hh"

namespace Amanzi {
namespace Operators {

class OperatorAccumulation {
 public:
  OperatorAccumulation(AmanziMesh::Entity_kind entity,
                       Teuchos::RCP<Operator> global_op) :
      global_op_(global_op),
      mesh_(Teuchos::null)
  {
    InitAccumulation_(entity);
  }

  OperatorAccumulation(AmanziMesh::Entity_kind entity,
                    Teuchos::RCP<AmanziMesh::Mesh> mesh) :
      global_op_(Teuchos::null),
      mesh_(mesh)
  {
    InitAccumulation_(entity);
  }

  OperatorAccumulation(AmanziMesh::Entity_kind entity,
                    Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      global_op_(Teuchos::null),
      mesh_(mesh)
  {
    InitAccumulation_(entity);
  }

  // update methods
  // -- update method for just adding to PC
  void AddAccumulationTerm(const Epetra_MultiVector& du);

  // -- linearized update methods with storage terms
  void AddAccumulationTerm(const CompositeVector& u0, const CompositeVector& s0, 
                           const CompositeVector& ss, double dT, const std::string& name);
  void AddAccumulationTerm(const CompositeVector& u0, const CompositeVector& ss, 
                           double dT, const std::string& name);
  void AddAccumulationTerm(const CompositeVector& u0, const CompositeVector& ss,
                           const std::string& name);

  // access (for developers only)
  int schema_dofs() { return local_op_schema_; }
  int schema_prec_dofs() { return global_op_schema_; }

  Teuchos::RCP<const Operator> global_operator() const { return global_op_; }
  Teuchos::RCP<Operator> global_operator() { return global_op_; }

 protected:
  void InitAccumulation_(AmanziMesh::Entity_kind entity);

 protected:
  // operator
  Teuchos::RCP<Operator> global_op_;
  Teuchos::RCP<Op> local_op_;
  int global_op_schema_, local_op_schema_;

  // mesh info
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int ncells_owned;
  int nfaces_owned;
  int nnodes_owned;

};


}  // namespace Operators
}  // namespace Amanzi


#endif


