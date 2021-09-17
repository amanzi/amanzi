/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "SaturationGasEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Simple constructor
****************************************************************** */
SaturationGasEvaluator::SaturationGasEvaluator(Teuchos::ParameterList& plist)
  : SecondaryVariableFieldEvaluator(plist)
{
  my_key_ = plist_.get<std::string>("my key");
  saturation_liquid_key_ = plist_.get<std::string>("saturation liquid key");
  dependencies_.insert(saturation_liquid_key_);
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> SaturationGasEvaluator::Clone() const {
  return Teuchos::rcp(new SaturationGasEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void SaturationGasEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& sl = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  auto& result_c = *result->ViewComponent("cell");

  int ncells = result_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = 1.0 - sl[0][c];
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void SaturationGasEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  auto& result_c = *result->ViewComponent("cell");

  int ncells = result_c.MyLength();
  if (wrt_key == saturation_liquid_key_) {
    for (int c = 0; c != ncells; ++c) result_c[0][c] = -1.0;
  }
}

}  // namespace Multiphase
}  // namespace Amanzi

