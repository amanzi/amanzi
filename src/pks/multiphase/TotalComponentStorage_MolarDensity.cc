/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Field evaluator for a total component storage (water, hydrogen,
  etc) storage, the conserved quantity:

    TCS = phi * (rho_l * s_l + rho_g * s_g)
*/

#include "TotalComponentStorage_MolarDensity.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Constructor.
****************************************************************** */
TotalComponentStorage_MolarDensity::TotalComponentStorage_MolarDensity(Teuchos::ParameterList& plist) :
    MultiphaseBaseEvaluator(plist) {
  Init_();
}


/* ******************************************************************
* Initialization.
****************************************************************** */
void TotalComponentStorage_MolarDensity::Init_()
{
  my_key_ = plist_.get<std::string>("my key");

  saturation_liquid_key_ = plist_.get<std::string>("saturation liquid key");
  porosity_key_ = plist_.get<std::string>("porosity key");
  molar_density_liquid_key_ = plist_.get<std::string>("molar density liquid key");
  molar_density_gas_key_ = plist_.get<std::string>("molar density gas key");

  dependencies_.insert(porosity_key_);
  dependencies_.insert(saturation_liquid_key_);
  dependencies_.insert(molar_density_liquid_key_);
  dependencies_.insert(molar_density_gas_key_);
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
TotalComponentStorage_MolarDensity::TotalComponentStorage_MolarDensity(
    const TotalComponentStorage_MolarDensity& other) : MultiphaseBaseEvaluator(other) {};


Teuchos::RCP<Evaluator> TotalComponentStorage_MolarDensity::Clone() const {
  return Teuchos::rcp(new TotalComponentStorage_MolarDensity(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void TotalComponentStorage_MolarDensity::Evaluate_(
    const State& S, const std::vector<CompositeVector*>& results)
{
  const auto& phi = *S->GetFieldData(porosity_key_)->ViewComponent("cell");
  const auto& sl = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  const auto& nl = *S->GetFieldData(molar_density_liquid_key_)->ViewComponent("cell");
  const auto& ng = *S->GetFieldData(molar_density_gas_key_)->ViewComponent("cell");

  auto& result_c = *result->ViewComponent("cell");
  int ncells = result->size("cell", false);

  for (int c = 0; c != ncells; ++c) {
    double tmpl = phi[0][c] * nl[0][c] * sl[0][c];
    double tmpg = phi[0][c] * ng[0][c] * (1.0 - sl[0][c]);
    result_c[0][c] = tmpl + tmpg;
  }
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void TotalComponentStorage_MolarDensity::EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& results)
{
  const auto& phi = *S.GetFieldData(porosity_key_)->ViewComponent("cell");
  const auto& sl = *S.GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  const auto& nl = *S.GetFieldData(molar_density_liquid_key_)->ViewComponent("cell");
  const auto& ng = *S.GetFieldData(molar_density_gas_key_)->ViewComponent("cell");

  auto& result_c = *result->ViewComponent("cell");
  int ncells = result->size("cell", false);

  if (wrt_key == porosity_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = sl[0][c] * nl[0][c] + (1.0 - sl[0][c]) * ng[0][c];
    }
  }
  else if (wrt_key == saturation_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * (nl[0][c] - ng[0][c]);
    }
  }

  else if (wrt_key == molar_density_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * sl[0][c];
    }
  } else if (wrt_key == molar_density_gas_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * (1.0 - sl[0][c]);
    }
  }
}

}  // namespace Multiphase
}  // namespace Amanzi
