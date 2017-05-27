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

#include "Epetra_MultiVector.h"

#include "BCsList.hh"
#include "PDE_Helper.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class Reaction : public BCsList, public PDE_Helper {
 public:
  Reaction(Teuchos::ParameterList& plist, Teuchos::RCP<Operator> global_op) :
      K_(Teuchos::null) {
    InitReaction_(plist);
  }

  Reaction(Teuchos::ParameterList& plist, Teuchos::RCP<const AmanziMesh::Mesh> mesh) : 
      K_(Teuchos::null),
      PDE_Helper(mesh) {
    InitReaction_(plist);
  }

  // required members 
  // -- setup
  virtual void Setup(Teuchos::RCP<Epetra_MultiVector>& K) { K_ = K; }
  // -- data
  virtual void UpdateMatrices(const CompositeVector& u);
  virtual void UpdateMatrices(const CompositeVector& u, const CompositeVector& dhdT) {};
  // -- results -- determine advected flux of u
  void UpdateFlux(const CompositeVector& h , const CompositeVector& u,
                  const Teuchos::RCP<BCs>& bc, CompositeVector& flux);
  
  // boundary conditions
  void ApplyBCs(bool primary, bool eliminate);

 private:
  void InitReaction_(Teuchos::ParameterList& plist);

 protected:
  Teuchos::RCP<const Epetra_MultiVector> K_;

 private:
  Schema global_schema_col_, global_schema_row_;
  Schema local_schema_col_, local_schema_row_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

