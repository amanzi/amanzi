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
    n_(0)
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
  const auto& s_l = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  const auto& n_l = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell");
  const auto& x_l = *S->GetFieldData("molar_fraction_liquid")->ViewComponent("cell");
  const auto& n_g = *S->GetFieldData("molar_density_gas")->ViewComponent("cell");

  auto& result_c = *result->ViewComponent("cell");
  int ncells = result->size("cell", false);

  for (int c = 0; c != ncells; ++c) {
    double tmpl = phi[0][c] * n_l[0][c] * s_l[0][c];
    double tmpg = phi[0][c] * n_g[0][c] * (1.0 - s_l[0][c]);
    result_c[0][c] = tmpl * x_l[n_][c] + tmpg * (1.0 - x_l[n_][c]);
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
  const auto& s_l = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  const auto& n_l = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell");
  const auto& x_l = *S->GetFieldData("molar_fraction_liquid")->ViewComponent("cell");
  const Epetra_MultiVector& n_g = *S->GetFieldData("molar_density_gas")->ViewComponent("cell");

  auto& result_c = *result->ViewComponent("cell");
  int ncells = result->size("cell", false);

  if (wrt_key == porosity_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = (s_l[0][c] * n_l[0][c] * x_l[n_][c] 
                     + (1.0 - s_l[0][c]) * n_g[0][c] * (1.0 - x_l[n_][c]));
    }
  }
  else if (wrt_key == saturation_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[n_][c] * (n_l[0][c] * x_l[n_][c] - n_g[0][c] * (1.0 - x_l[n_][c]));
    }
  }

  else if (wrt_key == "molar_density_liquid") {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * s_l[0][c] * x_l[n_][c];
    }
  } else if (wrt_key == "molar_density_gas") {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * (1.0 - s_l[0][c]) * (1.0 - x_l[n_][c]);
    }
  }

  else if (wrt_key == "molar_fraction_liquid") {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * (s_l[0][c] * n_l[0][c] - (1.0 - s_l[0][c]) * n_g[0][c]);
    }
  } else if (wrt_key == "molar_fraction_gas") {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * ((1.0 - s_l[0][c]) * n_g[0][c] - s_l[0][c] * n_l[0][c]);
    }
  } else {
    AMANZI_ASSERT(false);
  }
}

}  // namespace Multiphase
}  // namespace Amanzi
