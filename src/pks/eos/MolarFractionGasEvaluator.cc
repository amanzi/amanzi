/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Determining the molar fraction of a gas component within a gas mixture.
*/

#include "MolarFractionGasEvaluator.hh"
#include "SaturatedVaporPressureFactory.hh"

namespace Amanzi {
namespace AmanziEOS {

/* ******************************************************************
* Initialize various PK objects
****************************************************************** */
MolarFractionGasEvaluator::MolarFractionGasEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  // set up the actual model
  AMANZI_ASSERT(plist_.isSublist("vapor pressure model parameters"));
  SaturatedVaporPressureFactory svp_fac;
  svp_model_ = svp_fac.CreateVaporPressure(plist_.sublist("vapor pressure model parameters"));

  // process the list for my provided field.
  my_key_ = plist_.get<std::string>("molar fraction key");
  temp_key_ = plist_.get<std::string>("temperature key", "temperature");
  dependencies_.insert(temp_key_);
}


MolarFractionGasEvaluator::MolarFractionGasEvaluator(const MolarFractionGasEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    svp_model_(other.svp_model_),
    temp_key_(other.temp_key_) {};


Teuchos::RCP<FieldEvaluator> MolarFractionGasEvaluator::Clone() const {
  return Teuchos::rcp(new MolarFractionGasEvaluator(*this));
}


void MolarFractionGasEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  const double& p_atm = *(S->GetScalarData("atmospheric_pressure"));

  // evaluate p_s / p_atm
  for (CompositeVector::name_iterator comp=result->begin(); comp!=result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp, false));
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp, false));

    int count = result->size(*comp);
    for (int id = 0; id != count; ++id) {
      AMANZI_ASSERT(temp_v[0][id] > 200.0);
      result_v[0][id] = svp_model_->Pressure(temp_v[0][id]) / p_atm;
    }
  }
}


void MolarFractionGasEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S, Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result)
{
  AMANZI_ASSERT(wrt_key == temp_key_);

  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  const double& p_atm = *(S->GetScalarData("atmospheric_pressure"));

  for (auto comp = result->begin(); comp != result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp, false));
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp, false));

    int count = result->size(*comp);
    for (int id = 0; id != count; ++id) {
      result_v[0][id] = svp_model_->DPressureDT(temp_v[0][id]) / p_atm;
    }
  }
}

}  // namespace AmanziEOS
}  // namespace Amanzi

