/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Field evaluator for total water storage, the conserved quantity:

    TCS = phi * rho_l * s_l.
*/

#include "TotalWaterStorage.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Constructor.
****************************************************************** */
TotalWaterStorage::TotalWaterStorage(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
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

  dependencies_.insert(std::string(porosity_key_));
  dependencies_.insert(std::string(saturation_liquid_key_));
  dependencies_.insert(std::string("molar_density_liquid"));
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
TotalWaterStorage::TotalWaterStorage(const TotalWaterStorage& other) :
    SecondaryVariableFieldEvaluator(other) {};


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
  const auto& s_l = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  const auto& n_l = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell");

  auto& result_c = *result->ViewComponent("cell");
  int ncells = result->size("cell", false);

  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = phi[0][c] * s_l[0][c] * n_l[0][c];
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
  const auto& s_l = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  const auto& n_l = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell");

  auto& result_c = *result->ViewComponent("cell");
  int ncells = result->size("cell", false);

  if (wrt_key == porosity_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = s_l[0][c] * n_l[0][c];
    }
  } else if (wrt_key == saturation_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * n_l[0][c];
    }
  } else if (wrt_key == "molar_density_liquid") {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * s_l[0][c];
    }
  } else {
    AMANZI_ASSERT(0);
  }
}

}  // namespace Multiphase
}  // namespace Amanzi
