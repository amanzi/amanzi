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

#include "BCsList.hh"
#include "PDE_Helper.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class Advection : public BCsList, public PDE_Helper {
 public:
  Advection(Teuchos::ParameterList& plist,
            Teuchos::RCP<Operator> global_op) {
    mesh_ = Teuchos::null;
    global_op_ = global_op;
  }

  Advection(Teuchos::ParameterList& plist,
            Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      PDE_Helper(mesh) {
    global_op_ = Teuchos::null;
  }

  // main members
  // -- setup
  virtual void Setup(const CompositeVector& u) = 0;
  // -- populate data: arrays of local matrices
  virtual void UpdateMatrices(const CompositeVector& u) = 0;
  virtual void UpdateMatrices(const CompositeVector& u, const CompositeVector& dhdT) = 0;

  // -- results: determine advected flux of u
  virtual void UpdateFlux(const CompositeVector& h , const CompositeVector& u,
                          const Teuchos::RCP<BCs>& bc, CompositeVector& flux) = 0;
  
  // access
  Teuchos::RCP<Operator> global_operator() { return global_op_; }
  Teuchos::RCP<const Operator> global_operator() const { return global_op_; }

 protected:
  Schema global_schema_row_, global_schema_col_;
  Schema local_schema_col_, local_schema_row_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

