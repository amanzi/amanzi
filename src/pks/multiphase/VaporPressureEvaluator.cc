/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "errors.hh"

#include "EOSFactory.hh"
#include "EOS_SaturatedVaporPressure.hh"
#include "CommonDefs.hh"

#include "ModelMeshPartition.hh"
#include "MultiphaseDefs.hh"
#include "VaporPressureEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Simple constructor
****************************************************************** */
VaporPressureEvaluator::VaporPressureEvaluator(
    Teuchos::ParameterList& plist, Teuchos::RCP<WRMmpPartition> wrm)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist),
      wrm_(wrm)
{
  if (my_keys_.size() == 0) {
    my_keys_.push_back(std::make_pair(plist_.get<std::string>("my key"), Tags::DEFAULT));
  }
  temperature_key_ = plist_.get<std::string>("temperature key");
  molar_density_liquid_key_ = plist_.get<std::string>("molar density liquid key");
  saturation_liquid_key_ = plist_.get<std::string>("saturation liquid key");

  dependencies_.insert(std::make_pair(temperature_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(molar_density_liquid_key_, Tags::DEFAULT));
  dependencies_.insert(std::make_pair(saturation_liquid_key_, Tags::DEFAULT));

  AmanziEOS::EOSFactory<AmanziEOS::EOS_SaturatedVaporPressure> svp_factory;
  svp_ = svp_factory.Create(plist);
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<Evaluator> VaporPressureEvaluator::Clone() const {
  return Teuchos::rcp(new VaporPressureEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void VaporPressureEvaluator::Evaluate_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& temp = *S.Get<CompositeVector>(temperature_key_).ViewComponent("cell");
  const auto& sat = *S.Get<CompositeVector>(saturation_liquid_key_).ViewComponent("cell");
  const auto& nl = *S.Get<CompositeVector>(molar_density_liquid_key_).ViewComponent("cell");
  auto& result_c = *results[0]->ViewComponent("cell");

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
void VaporPressureEvaluator::EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results) 
{
  const auto& temp = *S.Get<CompositeVector>(temperature_key_).ViewComponent("cell");
  const auto& sat = *S.Get<CompositeVector>(saturation_liquid_key_).ViewComponent("cell");
  const auto& nl = *S.Get<CompositeVector>(molar_density_liquid_key_).ViewComponent("cell");
  auto& result_c = *results[0]->ViewComponent("cell");

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

