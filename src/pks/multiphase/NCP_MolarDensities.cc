/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Field evaluator for noninear complimentary problem, function G.
*/

#include "NCP_MolarDensities.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Constructor.
****************************************************************** */
NCP_MolarDensities::NCP_MolarDensities(Teuchos::ParameterList& plist)
  : MultiphaseBaseEvaluator(plist)
{
  my_key_ = plist_.get<std::string>("my key");
  x_vapor_key_ = plist_.get<std::string>("mole fraction vapor key");
  molar_density_gas_key_ = plist_.get<std::string>("molar density gas key");
  tcc_gas_key_ = plist_.get<std::string>("tcc gas key");

  dependencies_.insert(x_vapor_key_);
  dependencies_.insert(molar_density_gas_key_);
  dependencies_.insert(tcc_gas_key_);
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
NCP_MolarDensities::NCP_MolarDensities(const NCP_MolarDensities& other)
  : MultiphaseBaseEvaluator(other) {};


Teuchos::RCP<FieldEvaluator> NCP_MolarDensities::Clone() const {
  return Teuchos::rcp(new NCP_MolarDensities(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void NCP_MolarDensities::EvaluateField_(
    const Teuchos::Ptr<State>& S, const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& xg = *S->GetFieldData(x_vapor_key_)->ViewComponent("cell");
  const auto& ng = *S->GetFieldData(molar_density_gas_key_)->ViewComponent("cell");
  const auto& tcc = *S->GetFieldData(tcc_gas_key_)->ViewComponent("cell");

  auto& result_c = *result->ViewComponent("cell");
  int ncells = result->size("cell", false);

  for (int c = 0; c != ncells; ++c) {
    double sum = ng[0][c] * xg[0][c];
    for (int i = 0; i < tcc.NumVectors(); ++i) sum += tcc[i][c];
    result_c[0][c] = ng[0][c] - sum;
  }      
}


/* ******************************************************************
* Required member: field derivative calculation.
****************************************************************** */
void NCP_MolarDensities::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& xg = *S->GetFieldData(x_vapor_key_)->ViewComponent("cell");
  const auto& ng = *S->GetFieldData(molar_density_gas_key_)->ViewComponent("cell");

  auto& result_c = *result->ViewComponent("cell");
  int ncells = result->size("cell", false);

  if (wrt_key == x_vapor_key_) {
    for (int c = 0; c != ncells; ++c) result_c[0][c] = -ng[0][c];
  }
  else if (wrt_key == molar_density_gas_key_) {
    for (int c = 0; c != ncells; ++c) result_c[0][c] = 1.0 - xg[0][c];
  }
  else if (wrt_key == tcc_gas_key_) {
    for (int c = 0; c != ncells; ++c) result_c[0][c] = -1.0;
  }
}

}  // namespace Multiphase
}  // namespace Amanzi

