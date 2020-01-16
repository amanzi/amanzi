/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Field evaluator for a total component stirage (water, hydrogen,
  etc) storage, the conserved quantity:

    TCS = phi * (rho_l * s_l * X_l + rho_g * s_g * X_g)

  where x_p is the mole fraction of a component in phase p.
*/

#include "TotalComponentStorage.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Constructor.
****************************************************************** */
TotalComponentStorage::TotalComponentStorage(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist),
    n_(0),
    kH_(1.0e-2)
{
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

  dependencies_.insert(std::string(porosity_key_));
  dependencies_.insert(std::string(saturation_liquid_key_));
  dependencies_.insert(std::string("molar_density_liquid"));
  dependencies_.insert(std::string("molar_fraction_liquid"));
  dependencies_.insert(std::string("molar_density_gas"));
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
TotalComponentStorage::TotalComponentStorage(const TotalComponentStorage& other) :
    SecondaryVariableFieldEvaluator(other) {};


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
  const auto& nl = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell");
  const auto& xl = *S->GetFieldData("molar_fraction_liquid")->ViewComponent("cell");
  const auto& ng = *S->GetFieldData("molar_density_gas")->ViewComponent("cell");

  auto& result_c = *result->ViewComponent("cell");
  int ncells = result->size("cell", false);

  for (int c = 0; c != ncells; ++c) {
    double tmpl = phi[0][c] * nl[0][c] * sl[0][c];
    double tmpg = phi[0][c] * ng[0][c] * (1.0 - sl[0][c]);
    double xg = xl[n_][c] / kH_;
    result_c[0][c] = tmpl * xl[n_][c] + tmpg * xg;
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
  const auto& nl = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell");
  const auto& xl = *S->GetFieldData("molar_fraction_liquid")->ViewComponent("cell");
  const auto& ng = *S->GetFieldData("molar_density_gas")->ViewComponent("cell");

  auto& result_c = *result->ViewComponent("cell");
  int ncells = result->size("cell", false);

  if (wrt_key == porosity_key_) {
    for (int c = 0; c != ncells; ++c) {
      double xg = xl[n_][c] / kH_;
      result_c[0][c] = sl[0][c] * nl[0][c] * xl[n_][c] + (1.0 - sl[0][c]) * ng[0][c] * xg;
    }
  }
  else if (wrt_key == saturation_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      double xg = xl[n_][c] / kH_;
      result_c[0][c] = phi[n_][c] * (nl[0][c] * xl[n_][c] - ng[0][c] * xg);
    }
  }

  else if (wrt_key == "molar_density_liquid") {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * sl[0][c] * xl[n_][c];
    }
  } else if (wrt_key == "molar_density_gas") {
    for (int c = 0; c != ncells; ++c) {
      double xg = xl[n_][c] / kH_;
      result_c[0][c] = phi[0][c] * (1.0 - sl[0][c]) * xg;
    }
  }

  else if (wrt_key == "molar_fraction_liquid") {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * (sl[0][c] * nl[0][c] + (1.0 - sl[0][c]) * ng[0][c] / kH_);
    }
  }
}

}  // namespace Multiphase
}  // namespace Amanzi
