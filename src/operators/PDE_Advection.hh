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

#ifndef AMANZI_OPERATOR_PDE_ADVECTION_HH_
#define AMANZI_OPERATOR_PDE_ADVECTION_HH_

#include "PDE_HelperDiscretization.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "Schema.hh"


namespace Amanzi {
namespace Operators {

class PDE_Advection : public PDE_HelperDiscretization {
 public:
  PDE_Advection(Teuchos::ParameterList& plist,
                const Teuchos::RCP<Operator>& global_op) : PDE_HelperDiscretization(global_op) {};

  PDE_Advection(Teuchos::ParameterList& plist,
                const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : PDE_HelperDiscretization(mesh) {
    global_op_ = Teuchos::null;
  }

  virtual ~PDE_Advection() {};

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
                          const Teuchos::Ptr<CompositeVector>& u) = 0;
 protected:
  std::string name_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

