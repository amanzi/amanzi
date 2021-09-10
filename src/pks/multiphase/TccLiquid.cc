/*
  MultiPhase PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Calculates liquid mole fraction from partial gas pressure.
*/

#include "ModelMeshPartition.hh"
#include "MultiphaseDefs.hh"
#include "TccLiquid.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Simple constructor
****************************************************************** */
TccLiquid::TccLiquid(Teuchos::ParameterList& plist)
  : MultiphaseBaseEvaluator(plist)
{
  my_key_ = plist_.get<std::string>("my key");
  tcc_gas_key_ = plist_.get<std::string>("tcc gas key");
  dependencies_.insert(tcc_gas_key_);
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> TccLiquid::Clone() const {
  return Teuchos::rcp(new TccLiquid(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void TccLiquid::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& tcc = *S->GetFieldData(tcc_gas_key_)->ViewComponent("cell");
  auto& result_c = *result->ViewComponent("cell");
  int ncells = result_c.MyLength();

  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = tcc[0][c] * kH_;
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void TccLiquid::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  auto& result_c = *result->ViewComponent("cell");
  int ncells = result_c.MyLength();

  if (wrt_key == tcc_gas_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = kH_;
    }
  }
}

}  // namespace Multiphase
}  // namespace Amanzi

