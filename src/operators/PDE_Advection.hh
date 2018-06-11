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
               Teuchos::RCP<Operator> global_op) : PDE_HelperDiscretization(global_op) {};

  PDE_Advection(Teuchos::ParameterList& plist,
                Teuchos::RCP<const AmanziMesh::Mesh> mesh) : PDE_HelperDiscretization(mesh) {
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
                          const Teuchos::Ptr<CompositeVector>& u) = 0;
  
  // -- matrix modification
  //    primary=true indicates that the operator updates both matrix and right-hand
  //      side using BC data. If primary=false, only matrix is changed.
  //    eliminate=true indicates that we eliminate essential BCs for a trial 
  //      function, i.e. zeros go in the corresponding matrix columns and 
  //      right-hand side is modified using BC values. This is the optional 
  //      parameter that enforces symmetry for a symmetric tree  operators.
  //    essential_eqn=true indicates that the operator places a positive number on 
  //      the main matrix diagonal for the case of essential BCs. This is the
  //      implementtion trick.
  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn) override {};

 protected:
  Schema global_schema_row_, global_schema_col_;
  Schema local_schema_col_, local_schema_row_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

