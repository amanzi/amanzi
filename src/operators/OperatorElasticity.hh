/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_ELASTICITY_HH_
#define AMANZI_OPERATOR_ELASTICITY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "Tensor.hh"
#include "CompositeVector.hh"

#include "BCs.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class OperatorElasticity {
 public:
  OperatorElasticity(Teuchos::ParameterList& plist,
                     const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      plist_(plist),
      mesh_(mesh),
      K_(Teuchos::null),
      global_op_(Teuchos::null),
      ncells_owned(-1),
      ncells_wghost(-1),
      nfaces_owned(-1),
      nfaces_wghost(-1)
  {
    operator_type_ = OPERATOR_ELASTICITY;
    InitElasticity_(plist);
  }

  // main virtual members
  // -- setup 
  virtual void SetTensorCoefficient(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K);

  // -- creation of an operator
  virtual void UpdateMatrices();

  // -- matrix modification
  virtual void ApplyBCs(bool primary, bool eliminate);

  // boundary conditions (BC) require information on test and
  // trial spaces. For a single PDE, these BCs could be the same.
  virtual void SetBCs(const Teuchos::RCP<BCs>& bc_trial,
                      const Teuchos::RCP<BCs>& bc_test) {
    bcs_trial_.clear();
    bcs_test_.clear();

    bcs_trial_.push_back(bc_trial);
    bcs_test_.push_back(bc_test);
  }

  // access
  Teuchos::RCP<const Operator> global_operator() const { return global_op_; }
  Teuchos::RCP<Operator> global_operator() { return global_op_; }
  const Schema& global_schema_col() { return global_schema_col_; }
  const Schema& schema_col() { return local_schema_col_; }
  const Schema& schema_row() { return local_schema_row_; }

 protected:
  void InitElasticity_(Teuchos::ParameterList& plist);
  void ApplyBCs_Face_(const Teuchos::Ptr<BCs>& bc_f, bool primary, bool eliminate);
  void ApplyBCs_Node_(const Teuchos::Ptr<BCs>& bc_v, bool primary, bool eliminate);

 protected:
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K_;

  // operator
  Teuchos::RCP<Operator> global_op_;
  Teuchos::RCP<Op> local_op_;

  Schema global_schema_col_, global_schema_row_;
  Schema local_schema_col_, local_schema_row_;

  std::vector<Teuchos::RCP<BCs> > bcs_trial_, bcs_test_;
  OperatorType operator_type_;

  // mesh info
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;

  // miscaleneous
  Teuchos::ParameterList plist_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


