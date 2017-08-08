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

#include "PDE_Helper.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "Schema.hh"


namespace Amanzi {
namespace Operators {

class Advection : public PDE_Helper {
 public:
  Advection(Teuchos::ParameterList& plist,
            Teuchos::RCP<Operator> global_op) : PDE_Helper(global_op) {};

  Advection(Teuchos::ParameterList& plist,
            Teuchos::RCP<const AmanziMesh::Mesh> mesh) : PDE_Helper(mesh) {
    global_op_ = Teuchos::null;
  }

  // main members
  // -- setup
  virtual void Setup(const CompositeVector& u) = 0;

  // -- standard interface for flux calculation
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& p,
                          const Teuchos::Ptr<CompositeVector>& u) override {};

  // -- extended interface for flux calculation
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& h,
                          const Teuchos::Ptr<const CompositeVector>& p,
                          const Teuchos::RCP<BCs>& bc, 
                          Teuchos::Ptr<CompositeVector>& u) = 0;
  
 protected:
  Schema global_schema_row_, global_schema_col_;
  Schema local_schema_col_, local_schema_row_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

