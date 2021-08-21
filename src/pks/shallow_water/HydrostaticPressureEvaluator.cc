/*
  Shallow Water PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/

#include "errors.hh"

#include "HydrostaticPressureEvaluator.hh"

namespace Amanzi {
namespace ShallowWater {

/* ******************************************************************
* Simple constructor
****************************************************************** */
HydrostaticPressureEvaluator::HydrostaticPressureEvaluator(Teuchos::ParameterList& plist)
  : SecondaryVariableFieldEvaluator(plist)
{
  my_key_ = "surface-ponded_pressure";
  ponded_depth_key_ = plist_.get<std::string>("ponded depth key", "surface-ponded_depth");

  dependencies_.insert(ponded_depth_key_);
}


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> HydrostaticPressureEvaluator::Clone() const {
  return Teuchos::rcp(new HydrostaticPressureEvaluator(*this));
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void HydrostaticPressureEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& h_c = *S->GetFieldData(ponded_depth_key_)->ViewComponent("cell");
  const double rho = *S->GetScalarData("const_fluid_density");
  const double patm = *S->GetScalarData("atmospheric_pressure");
  
  double tmp[1];
  S->GetConstantVectorData("gravity", "state")->Norm2(tmp);
  double g = tmp[0];
  
  auto& result_c = *result->ViewComponent("cell");

  int ncells = result_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    result_c[0][c] = rho * g * h_c[0][c] + patm;
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void HydrostaticPressureEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  const auto& h_c = *S->GetFieldData(ponded_depth_key_)->ViewComponent("cell");
  const double rho = *S->GetScalarData("const_fluid_density");
  
  double tmp[1];
  S->GetConstantVectorData("gravity", "state")->Norm2(tmp);
  double g = tmp[0];
  
  auto& result_c = *result->ViewComponent("cell");

  int ncells = result_c.MyLength();
  if (wrt_key == ponded_depth_key_) {
    for (int c = 0; c != ncells; ++c) {
      result_c[0][c] = rho * g;
    }
  } else {
    AMANZI_ASSERT(false);
  }
}

}  // namespace ShallowWater
}  // namespace Amanzi

