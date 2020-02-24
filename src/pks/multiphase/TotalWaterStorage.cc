/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Field evaluator for total water storage, the conserved quantity:

    TCS = phi * (eta_l * s_l + eta_g * s_g * x_v)
*/

#include "TotalWaterStorage.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Constructor.
****************************************************************** */
TotalWaterStorage::TotalWaterStorage(Teuchos::ParameterList& plist) :
    MultiphaseBaseEvaluator(plist)
{
  Init_();
}


/* ******************************************************************
* Initialization.
****************************************************************** */
void TotalWaterStorage::Init_()
{
  my_key_ = plist_.get<std::string>("my key");
  saturation_liquid_key_ = plist_.get<std::string>("saturation liquid key");
  porosity_key_ = plist_.get<std::string>("porosity key");
  molar_density_liquid_key_ = plist_.get<std::string>("molar density liquid key");
  molar_density_gas_key_ = plist_.get<std::string>("molar density gas key");
  x_vapor_key_ = plist_.get<std::string>("mole fraction vapor key");

  dependencies_.insert(porosity_key_);
  dependencies_.insert(saturation_liquid_key_);
  dependencies_.insert(molar_density_liquid_key_);
  dependencies_.insert(molar_density_gas_key_);
  dependencies_.insert(x_vapor_key_);
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
TotalWaterStorage::TotalWaterStorage(const TotalWaterStorage& other) :
    MultiphaseBaseEvaluator(other) {};


Teuchos::RCP<FieldEvaluator> TotalWaterStorage::Clone() const {
  return Teuchos::rcp(new TotalWaterStorage(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void TotalWaterStorage::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& phi = *S->GetFieldData(porosity_key_)->ViewComponent("cell");
  const auto& sl = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  const auto& nl = *S->GetFieldData(molar_density_liquid_key_)->ViewComponent("cell");
  const auto& ng = *S->GetFieldData(molar_density_gas_key_)->ViewComponent("cell");
  const auto& vg = *S->GetFieldData(x_vapor_key_)->ViewComponent("cell");

  auto& result_c = *result->ViewComponent("cell");
  int ncells = result->size("cell", false);

  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = phi[0][c] * (sl[0][c] * nl[0][c]
                   + (1.0 - sl[0][c]) * ng[0][c] * vg[0][c]);
  }      
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void TotalWaterStorage::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& phi = *S->GetFieldData(porosity_key_)->ViewComponent("cell");
  const auto& sl = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  const auto& nl = *S->GetFieldData(molar_density_liquid_key_)->ViewComponent("cell");
  const auto& ng = *S->GetFieldData(molar_density_gas_key_)->ViewComponent("cell");
  const auto& vg = *S->GetFieldData(x_vapor_key_)->ViewComponent("cell");

  auto& result_c = *result->ViewComponent("cell");
  int ncells = result->size("cell", false);

  if (wrt_key == porosity_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = sl[0][c] * nl[0][c] + (1.0 - sl[0][c]) * ng[0][c] * vg[0][c];
    }
  } else if (wrt_key == saturation_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * (nl[0][c] - ng[0][c] * vg[0][c]);
    }
  } else if (wrt_key == molar_density_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * sl[0][c];
    }
  } else if (wrt_key == molar_density_gas_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * (1.0 - sl[0][c]) * vg[0][c];
    }
  } else if (wrt_key == x_vapor_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * (1.0 - sl[0][c]) * ng[0][c];
    }
  }
}

}  // namespace Multiphase
}  // namespace Amanzi
