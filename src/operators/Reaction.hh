/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Reaction operator. 
*/

#ifndef AMANZI_OPERATOR_REACTION_HH_
#define AMANZI_OPERATOR_REACTION_HH_

#include <string>

#include "BCsList.hh"
#include "PDE_Helper.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class Reaction : public BCsList, public PDE_Helper {
 public:
  Reaction(Teuchos::ParameterList& plist, Teuchos::RCP<Operator> global_op) {
    InitReaction_(plist);
  }

  Reaction(Teuchos::ParameterList& plist,
           Teuchos::RCP<const AmanziMesh::Mesh> mesh) : PDE_Helper(mesh) {
    InitReaction_(plist);
  }

  // required members 
  // -- setup
  virtual void Setup(const CompositeVector& u) {};
  // -- data
  virtual void UpdateMatrices(const CompositeVector& u);
  virtual void UpdateMatrices(const CompositeVector& u, const CompositeVector& dhdT) {};
  // -- results -- determine advected flux of u
  void UpdateFlux(const CompositeVector& h , const CompositeVector& u,
                  const Teuchos::RCP<BCs>& bc, CompositeVector& flux);
  
  // boundary conditions
  void ApplyBCs(bool primary, bool eliminate);

  // access
  Teuchos::RCP<const Operator> global_operator() const { return global_op_; }
  Teuchos::RCP<Operator> global_operator() { return global_op_; }

 private:
  void InitReaction_(Teuchos::ParameterList& plist);

 private:
  Schema global_schema_col_, global_schema_row_;
  Schema local_schema_col_, local_schema_row_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

