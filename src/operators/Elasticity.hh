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

#include "BCsList.hh"
#include "PDE_Helper.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"
#include "Schema.hh"

namespace Amanzi {
namespace Operators {

class Elasticity : public BCsList, public PDE_Helper {
 public:
  Elasticity(Teuchos::ParameterList& plist,
             const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      PDE_Helper(mesh),
      plist_(plist),
      K_(Teuchos::null),
      K_default_(1.0)
  {
    global_op_ = Teuchos::null;
    operator_type_ = OPERATOR_ELASTICITY;
    InitElasticity_(plist);
  }

  // main virtual members
  // -- setup 
  void SetTensorCoefficient(const Teuchos::RCP<std::vector<WhetStone::Tensor> >& K);
  void SetTensorCoefficient(double K);

  // -- creation of an operator
  virtual void UpdateMatrices();

  // -- matrix modification
  virtual void ApplyBCs(bool primary, bool eliminate);

  // access
  Teuchos::RCP<const Operator> global_operator() const { return global_op_; }
  Teuchos::RCP<Operator> global_operator() { return global_op_; }
  const Schema& global_schema_col() { return global_schema_col_; }
  const Schema& schema_col() { return local_schema_col_; }
  const Schema& schema_row() { return local_schema_row_; }

 protected:
  void InitElasticity_(Teuchos::ParameterList& plist);

 protected:
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K_;
  double K_default_;

  // operator and schemas
  Schema global_schema_col_, global_schema_row_;
  Schema local_schema_col_, local_schema_row_;

  OperatorType operator_type_;

  // miscaleneous
  Teuchos::ParameterList plist_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif


