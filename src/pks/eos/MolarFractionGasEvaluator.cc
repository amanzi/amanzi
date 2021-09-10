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

#include "EOSFactory.hh"

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
  EOSFactory<EOS_SaturatedVaporPressure> svp_fac;
  svp_model_ = svp_fac.CreateEOS(plist_.sublist("vapor pressure model parameters"));

  // process the list for my provided field.
  if (my_key_ == "")
    my_key_ = plist_.get<std::string>("molar fraction key");

  // set up dependencies
  std::string domain = Keys::getDomain(my_key_);
  temp_key_ = plist_.get<std::string>("temperature key", Keys::getKey(domain, "temperature"));
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
  double p_atm = *(S->GetScalarData("atmospheric_pressure"));

  // evaluate p_s / p_atm
  for (auto comp = result->begin(); comp != result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp));
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp));

    int count = result->size(*comp);
    for (int i = 0; i != count; ++i) {
      result_v[0][i] = svp_model_->Pressure(temp_v[0][i]) / p_atm;
    }
  }
}


void MolarFractionGasEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S, Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result)
{
  AMANZI_ASSERT(wrt_key == temp_key_);

  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  double p_atm = *(S->GetScalarData("atmospheric_pressure"));

  for (auto comp = result->begin(); comp != result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp));
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp));

    int count = result->size(*comp);
    for (int i = 0; i != count; ++i) {
      result_v[0][i] = svp_model_->DPressureDT(temp_v[0][i]) / p_atm;
    }
  }
}

}  // namespace AmanziEOS
}  // namespace Amanzi

