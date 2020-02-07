/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "ModelMeshPartition.hh"
#include "MultiphaseDefs.hh"
#include "MoleFractionLiquid.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Simple constructor
****************************************************************** */
MoleFractionLiquid::MoleFractionLiquid(Teuchos::ParameterList& plist)
  : MultiphaseBaseEvaluator(plist)
{
  my_key_ = plist_.get<std::string>("my key");
  pressure_gas_key_ = plist_.get<std::string>("pressure gas key");
  x_gas_key_ = plist_.get<std::string>("mole fraction gas key");

  dependencies_.insert(std::string(pressure_gas_key_));
  dependencies_.insert(std::string(x_gas_key_));
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> MoleFractionLiquid::Clone() const {
  return Teuchos::rcp(new MoleFractionLiquid(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void MoleFractionLiquid::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& pg_c = *S->GetFieldData(pressure_gas_key_)->ViewComponent("cell");
  const auto& xg_c = *S->GetFieldData(x_gas_key_)->ViewComponent("cell");
  auto& result_c = *result->ViewComponent("cell");

  int ncells = result_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = pg_c[0][c] * xg_c[n_][c] * kH_;
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void MoleFractionLiquid::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& pg_c = *S->GetFieldData(pressure_gas_key_)->ViewComponent("cell");
  const auto& xg_c = *S->GetFieldData(x_gas_key_)->ViewComponent("cell");
  auto& result_c = *result->ViewComponent("cell");
 
  int ncells = result_c.MyLength();
  if (wrt_key == pressure_gas_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = xg_c[n_][c] * kH_;
    }
  } else if (wrt_key == x_gas_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = pg_c[0][c] * kH_;
    }
  }
}

}  // namespace Multiphase
}  // namespace Amanzi

