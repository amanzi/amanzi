/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Field evaluator for noninear complimentary problem, function G.
*/

#include "NCP_MoleFractions.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Constructor.
****************************************************************** */
NCP_MoleFractions::NCP_MoleFractions(Teuchos::ParameterList& plist)
  : MultiphaseBaseEvaluator(plist)
{
  my_key_ = plist_.get<std::string>("my key");
  x_vapor_key_ = plist_.get<std::string>("mole fraction vapor key");
  x_gas_key_ = plist_.get<std::string>("mole fraction gas key");

  dependencies_.insert(x_vapor_key_);
  dependencies_.insert(x_gas_key_);
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
NCP_MoleFractions::NCP_MoleFractions(const NCP_MoleFractions& other)
  : MultiphaseBaseEvaluator(other) {};


Teuchos::RCP<FieldEvaluator> NCP_MoleFractions::Clone() const {
  return Teuchos::rcp(new NCP_MoleFractions(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void NCP_MoleFractions::EvaluateField_(
    const Teuchos::Ptr<State>& S, const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& vg = *S->GetFieldData(x_vapor_key_)->ViewComponent("cell");
  const auto& xg = *S->GetFieldData(x_gas_key_)->ViewComponent("cell");

  auto& result_c = *result->ViewComponent("cell");
  int ncells = result->size("cell", false);

  for (int c = 0; c != ncells; ++c) {
    double sum(vg[0][c]);
    for (int i = 0; i < xg.NumVectors(); ++i) sum += xg[i][c];
    result_c[0][c] = 1.0 - sum;
  }      
}


/* ******************************************************************
* Required member: field derivative calculation.
****************************************************************** */
void NCP_MoleFractions::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  auto& result_c = *result->ViewComponent("cell");
  int ncells = result->size("cell", false);

  if (wrt_key == x_vapor_key_) {
    for (int c = 0; c != ncells; ++c) result_c[0][c] = -1.0;
  }
  else if (wrt_key == x_gas_key_) {
    for (int c = 0; c != ncells; ++c) result_c[0][c] = -1.0;
  }
}

}  // namespace Multiphase
}  // namespace Amanzi

