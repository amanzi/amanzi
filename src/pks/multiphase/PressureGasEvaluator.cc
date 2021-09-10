/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "errors.hh"

#include "ModelMeshPartition.hh"
#include "MultiphaseDefs.hh"
#include "PressureGasEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Simple constructor
****************************************************************** */
PressureGasEvaluator::PressureGasEvaluator(
    Teuchos::ParameterList& plist, Teuchos::RCP<WRMmpPartition> wrm)
  : SecondaryVariableFieldEvaluator(plist),
    wrm_(wrm)
{
  my_key_ = plist_.get<std::string>("my key");
  pressure_liquid_key_ = plist_.get<std::string>("pressure liquid key");
  saturation_liquid_key_ = plist_.get<std::string>("saturation liquid key");

  dependencies_.insert(std::string(pressure_liquid_key_));
  dependencies_.insert(std::string(saturation_liquid_key_));
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> PressureGasEvaluator::Clone() const {
  return Teuchos::rcp(new PressureGasEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void PressureGasEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& p_c = *S->GetFieldData(pressure_liquid_key_)->ViewComponent("cell");
  const auto& sat_c = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  auto& result_c = *result->ViewComponent("cell");

  int ncells = result_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = p_c[0][c] + wrm_->second[(*wrm_->first)[c]]->capillaryPressure(sat_c[0][c]);
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void PressureGasEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& sat_c = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  auto& result_c = *result->ViewComponent("cell");

  int ncells = result_c.MyLength();
  if (wrt_key == pressure_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = 1.0;
    }
  } else if (wrt_key == saturation_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = wrm_->second[(*wrm_->first)[c]]->dPc_dS(sat_c[0][c]);
    }
  }
}

}  // namespace Multiphase
}  // namespace Amanzi

