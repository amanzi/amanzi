/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Molar mobility is defined as the product of molar density and mobility.
*/

#include "MultiphaseTypeDefs.hh"
#include "MolarMobilityEvaluator.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Two constructors.
****************************************************************** */
MolarMobilityEvaluator::MolarMobilityEvaluator(Teuchos::ParameterList& plist,
                                               const Teuchos::RCP<WRMmpPartition>& wrm)
  : SecondaryVariableFieldEvaluator(plist),
    wrm_(wrm)
{
  phase_name_ = plist.get<std::string>("phase name");
  my_key_ = plist.get<std::string>("my key");
  saturation_liquid_key_ = plist.get<std::string>("saturation key", "saturation_liquid");
  viscosity_key_ = plist.get<std::string>("viscosity key");
  molar_density_key_ = plist.get<std::string>("molar density key");

  dependencies_.insert(saturation_liquid_key_);
  dependencies_.insert(viscosity_key_);
}


MolarMobilityEvaluator::MolarMobilityEvaluator(const MolarMobilityEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    wrm_(other.wrm_) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> MolarMobilityEvaluator::Clone() const {
  return Teuchos::rcp(new MolarMobilityEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void MolarMobilityEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& sat_c = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  const auto& vis_c = *S->GetFieldData(viscosity_key_)->ViewComponent("cell");
  const auto& eta_c = *S->GetFieldData(molar_density_key_)->ViewComponent("cell");
  auto& result_c = *result->ViewComponent("cell");

  int ncells = result_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    double kr = wrm_->second[(*wrm_->first)[c]]->k_relative(sat_c[0][c], phase_name_);
    result_c[0][c] = kr * eta_c[0][c] / vis_c[0][c];
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void MolarMobilityEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& sat_c = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  const auto& vis_c = *S->GetFieldData(viscosity_key_)->ViewComponent("cell");
  const auto& eta_c = *S->GetFieldData(molar_density_key_)->ViewComponent("cell");
  auto& result_c = *result->ViewComponent("cell");

  int ncells = result_c.MyLength();
  if (wrt_key == saturation_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      double dkr_ds = wrm_->second[(*wrm_->first)[c]]->dKdS(sat_c[0][c], phase_name_);
      result_c[0][c] = dkr_ds * eta_c[0][c] / vis_c[0][c];
    }
  } else if (wrt_key == viscosity_key_) {
    for (int c = 0; c != ncells; ++c) {
      double kr = wrm_->second[(*wrm_->first)[c]]->k_relative(sat_c[0][c], phase_name_);
      result_c[0][c] = -kr * eta_c[0][c] / (vis_c[0][c] * vis_c[0][c]);
    }
  } else if (wrt_key == molar_density_key_) {
    for (int c = 0; c != ncells; ++c) {
      double kr = wrm_->second[(*wrm_->first)[c]]->k_relative(sat_c[0][c], phase_name_);
      result_c[0][c] = kr / vis_c[0][c];
    }
  }
}

}  // namespace Flow
}  // namespace Amanzi
