/*
  Shallow Water PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "errors.hh"

#include "DischargeEvaluator.hh"

namespace Amanzi {
namespace ShallowWater {

/* ******************************************************************
* Simple constructor
****************************************************************** */
DischargeEvaluator::DischargeEvaluator(Teuchos::ParameterList& plist)
  : SecondaryVariableFieldEvaluator(plist)
{
  my_key_ = plist_.get<std::string>("my key", "surface-discharge");
  std::string domain = Keys::getDomain(my_key_);

  ponded_depth_key_ = plist_.get<std::string>("ponded depth key", Keys::getKey(domain, "ponded_depth"));
  velocity_key_ = plist_.get<std::string>("velocity key", Keys::getKey(domain, "velocity"));

  dependencies_.insert(std::string(ponded_depth_key_));
  dependencies_.insert(std::string(velocity_key_));
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> DischargeEvaluator::Clone() const {
  return Teuchos::rcp(new DischargeEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void DischargeEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& h_c = *S->GetFieldData(ponded_depth_key_)->ViewComponent("cell");
  const auto& u_c = *S->GetFieldData(velocity_key_)->ViewComponent("cell");
  auto& result_c = *result->ViewComponent("cell");

  int ncells = result_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    for (int i = 0; i < 2; ++i) {
      result_c[i][c] = h_c[0][c] * u_c[i][c];
    }
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void DischargeEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& h_c = *S->GetFieldData(ponded_depth_key_)->ViewComponent("cell");
  const auto& u_c = *S->GetFieldData(velocity_key_)->ViewComponent("cell");
  auto& result_c = *result->ViewComponent("cell");

  int ncells = result_c.MyLength();
  if (wrt_key == ponded_depth_key_) {
    for (int c = 0; c != ncells; ++c) {
      for (int i = 0; i < 2; ++i) {
        result_c[i][c] = u_c[i][c];
      }
    }
  } else {
    AMANZI_ASSERT(false);
  }
}

}  // namespace ShallowWater
}  // namespace Amanzi

