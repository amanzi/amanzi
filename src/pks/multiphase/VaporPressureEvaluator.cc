/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "errors.hh"

#include "CommonDefs.hh"

#include "ModelMeshPartition.hh"
#include "MultiphaseDefs.hh"
#include "SaturatedVaporPressureFactory.hh"
#include "VaporPressureEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Simple constructor
****************************************************************** */
VaporPressureEvaluator::VaporPressureEvaluator(
    Teuchos::ParameterList& plist, Teuchos::RCP<WRMmpPartition> wrm)
  : SecondaryVariableFieldEvaluator(plist),
    wrm_(wrm)
{
  my_key_ = plist_.get<std::string>("my key");
  temperature_key_ = plist_.get<std::string>("temperature key");
  molar_density_liquid_key_ = plist_.get<std::string>("molar density liquid key");
  saturation_liquid_key_ = plist_.get<std::string>("saturation liquid key");

  dependencies_.insert(temperature_key_);
  dependencies_.insert(molar_density_liquid_key_);
  dependencies_.insert(saturation_liquid_key_);

  AmanziEOS::SaturatedVaporPressureFactory svp_factory;
  svp_ = svp_factory.CreateVaporPressure(plist);
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> VaporPressureEvaluator::Clone() const {
  return Teuchos::rcp(new VaporPressureEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void VaporPressureEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& temp = *S->GetFieldData(temperature_key_)->ViewComponent("cell");
  const auto& sat = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  const auto& nl = *S->GetFieldData(molar_density_liquid_key_)->ViewComponent("cell");
  auto& result_c = *result->ViewComponent("cell");

  double R = CommonDefs::IDEAL_GAS_CONSTANT_R;

  int ncells = result_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    double pc = wrm_->second[(*wrm_->first)[c]]->capillaryPressure(sat[0][c]);
    result_c[0][c] = svp_->Pressure(temp[0][c]) * std::exp(-pc / (nl[0][c] * R * temp[0][c]));
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void VaporPressureEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& temp = *S->GetFieldData(temperature_key_)->ViewComponent("cell");
  const auto& sat = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  const auto& nl = *S->GetFieldData(molar_density_liquid_key_)->ViewComponent("cell");
  auto& result_c = *result->ViewComponent("cell");

  double R = CommonDefs::IDEAL_GAS_CONSTANT_R;

  int ncells = result_c.MyLength();
  if (wrt_key == temperature_key_) {
    for (int c = 0; c != ncells; ++c) {
      double pc = wrm_->second[(*wrm_->first)[c]]->capillaryPressure(sat[0][c]);
      double a = nl[0][c] * R * temp[0][c];
      double tmp = std::exp(-pc / a);
      result_c[0][c] = tmp * (svp_->DPressureDT(temp[0][c])
                            + svp_->Pressure(temp[0][c]) * pc / (a * temp[0][c]));
    }
  } else if (wrt_key == molar_density_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      double pc = wrm_->second[(*wrm_->first)[c]]->capillaryPressure(sat[0][c]);
      double a = nl[0][c] * R * temp[0][c];
      double tmp = std::exp(-pc / a);
      result_c[0][c] = svp_->Pressure(temp[0][c]) * tmp * pc / (a * nl[0][c]);
    }
  } else if (wrt_key == saturation_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      double pc = wrm_->second[(*wrm_->first)[c]]->capillaryPressure(sat[0][c]);
      double dPc_dS = wrm_->second[(*wrm_->first)[c]]->dPc_dS(sat[0][c]);
      double a = nl[0][c] * R * temp[0][c];
      double tmp = std::exp(-pc / a);
      result_c[0][c] = -svp_->Pressure(temp[0][c]) * tmp * dPc_dS / a;
    }
  }
}

}  // namespace Multiphase
}  // namespace Amanzi

