/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Field evaluator for a total component stirage (water, hydrogen,
  etc) storage, the conserved quantity:

    TCS = phi * (eta_l * s_l * x_l + eta_g * s_g * x_g)

  where X_p is the mole fraction of a component in phase p.
*/

#include "TotalComponentStorage.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Constructor.
****************************************************************** */
TotalComponentStorage::TotalComponentStorage(Teuchos::ParameterList& plist) :
    MultiphaseBaseEvaluator(plist) {
  Init_();
}


/* ******************************************************************
* Initialization.
****************************************************************** */
void TotalComponentStorage::Init_()
{
  my_key_ = plist_.get<std::string>("my key");

  saturation_liquid_key_ = plist_.get<std::string>("saturation liquid key");
  porosity_key_ = plist_.get<std::string>("porosity key");
  molar_density_liquid_key_ = plist_.get<std::string>("molar density liquid key");
  molar_density_gas_key_ = plist_.get<std::string>("molar density gas key");
  mole_fraction_liquid_key_ = plist_.get<std::string>("mole fraction liquid key");
  mole_fraction_gas_key_ = plist_.get<std::string>("mole fraction gas key");

  dependencies_.insert(porosity_key_);
  dependencies_.insert(saturation_liquid_key_);
  dependencies_.insert(molar_density_liquid_key_);
  dependencies_.insert(molar_density_gas_key_);
  dependencies_.insert(mole_fraction_liquid_key_);
  dependencies_.insert(mole_fraction_gas_key_);
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
TotalComponentStorage::TotalComponentStorage(const TotalComponentStorage& other) :
    MultiphaseBaseEvaluator(other) {};


Teuchos::RCP<FieldEvaluator> TotalComponentStorage::Clone() const {
  return Teuchos::rcp(new TotalComponentStorage(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void TotalComponentStorage::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& phi = *S->GetFieldData(porosity_key_)->ViewComponent("cell");
  const auto& sl = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  const auto& nl = *S->GetFieldData(molar_density_liquid_key_)->ViewComponent("cell");
  const auto& ng = *S->GetFieldData(molar_density_gas_key_)->ViewComponent("cell");
  const auto& xl = *S->GetFieldData(mole_fraction_liquid_key_)->ViewComponent("cell");
  const auto& xg = *S->GetFieldData(mole_fraction_gas_key_)->ViewComponent("cell");

  auto& result_c = *result->ViewComponent("cell");
  int ncells = result->size("cell", false);

  for (int c = 0; c != ncells; ++c) {
    double tmpl = phi[0][c] * nl[0][c] * sl[0][c];
    double tmpg = phi[0][c] * ng[0][c] * (1.0 - sl[0][c]);
    result_c[0][c] = tmpl * xl[0][c] + tmpg * xg[n_][c];
  }
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void TotalComponentStorage::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& phi = *S->GetFieldData(porosity_key_)->ViewComponent("cell");
  const auto& sl = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  const auto& nl = *S->GetFieldData(molar_density_liquid_key_)->ViewComponent("cell");
  const auto& ng = *S->GetFieldData(molar_density_gas_key_)->ViewComponent("cell");
  const auto& xl = *S->GetFieldData(mole_fraction_liquid_key_)->ViewComponent("cell");
  const auto& xg = *S->GetFieldData(mole_fraction_gas_key_)->ViewComponent("cell");

  auto& result_c = *result->ViewComponent("cell");
  int ncells = result->size("cell", false);

  if (wrt_key == porosity_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = sl[0][c] * nl[0][c] * xl[0][c] + (1.0 - sl[0][c]) * ng[0][c] * xg[n_][c];
    }
  }
  else if (wrt_key == saturation_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * (nl[0][c] * xl[0][c] - ng[0][c] * xg[n_][c]);
    }
  }

  else if (wrt_key == molar_density_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * sl[0][c] * xl[0][c];
    }
  } else if (wrt_key == molar_density_gas_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * (1.0 - sl[0][c]) * xg[n_][c];
    }
  }

  else if (wrt_key == mole_fraction_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * sl[0][c] * nl[0][c];
    }
  }
  else if (wrt_key == mole_fraction_gas_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * (1.0 - sl[0][c]) * ng[0][c];
    }
  }
}

}  // namespace Multiphase
}  // namespace Amanzi
