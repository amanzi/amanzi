/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Field evaluator for a total component storage (water, hydrogen,
  etc) storage, the conserved quantity:

    TCS = phi * (rho_l * s_l  rho_g * s_g)
*/

#include "TotalComponentStorageTest.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Constructor.
****************************************************************** */
TotalComponentStorageTest::TotalComponentStorageTest(Teuchos::ParameterList& plist) :
    MultiphaseBaseEvaluator(plist) {
  Init_();
}


/* ******************************************************************
* Initialization.
****************************************************************** */
void TotalComponentStorageTest::Init_()
{
  my_key_ = plist_.get<std::string>("my key");

  saturation_liquid_key_ = plist_.get<std::string>("saturation liquid key");
  porosity_key_ = plist_.get<std::string>("porosity key");

  dependencies_.insert(std::string(porosity_key_));
  dependencies_.insert(std::string(saturation_liquid_key_));
  dependencies_.insert(std::string("molar_density_liquid"));
  dependencies_.insert(std::string("molar_density_gas"));
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
TotalComponentStorageTest::TotalComponentStorageTest(const TotalComponentStorageTest& other)
  : MultiphaseBaseEvaluator(other) {};


Teuchos::RCP<FieldEvaluator> TotalComponentStorageTest::Clone() const {
  return Teuchos::rcp(new TotalComponentStorageTest(*this));
}


/* ******************************************************************
* Required member: field calculation.
****************************************************************** */
void TotalComponentStorageTest::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& phi = *S->GetFieldData(porosity_key_)->ViewComponent("cell");
  const auto& sl = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  const auto& nl = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell");
  const auto& ng = *S->GetFieldData("molar_density_gas")->ViewComponent("cell");

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
void TotalComponentStorageTest::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& phi = *S->GetFieldData(porosity_key_)->ViewComponent("cell");
  const auto& sl = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell");
  const auto& nl = *S->GetFieldData("molar_density_liquid")->ViewComponent("cell");
  const auto& ng = *S->GetFieldData("molar_density_gas")->ViewComponent("cell");

  auto& result_c = *result->ViewComponent("cell");
  int ncells = result->size("cell", false);

  if (wrt_key == porosity_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = sl[0][c] * nl[0][c] + (1.0 - sl[0][c]) * ng[0][c];
    }
  }
  else if (wrt_key == saturation_liquid_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[n_][c] * (nl[0][c] - ng[0][c]);
    }
  }

  else if (wrt_key == "molar_density_liquid") {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * sl[0][c];
    }
  } else if (wrt_key == "molar_density_gas") {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = phi[0][c] * (1.0 - sl[0][c]);
    }
  }
}

}  // namespace Multiphase
}  // namespace Amanzi
