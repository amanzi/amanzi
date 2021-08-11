/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  The WRM Evaluator simply calls the WRM with the correct arguments.
*/

#include "FlowDefs.hh"
#include "WRMEvaluator.hh"
#include "WRMPartition.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Two constructors.
****************************************************************** */
WRMEvaluator::WRMEvaluator(Teuchos::ParameterList& plist,
                           const Teuchos::RCP<WRMPartition>& wrm) :
    SecondaryVariableFieldEvaluator(plist),
    wrm_(wrm)
{
  InitializeFromPlist_();
}


WRMEvaluator::WRMEvaluator(const WRMEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    pressure_key_(other.pressure_key_),
    wrm_(other.wrm_) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> WRMEvaluator::Clone() const {
  return Teuchos::rcp(new WRMEvaluator(*this));
}


/* ******************************************************************
* Initialization.
****************************************************************** */
void WRMEvaluator::InitializeFromPlist_() {
  my_key_ = plist_.get<std::string>("saturation key", "saturation_liquid");

  // my dependency is pressure.
  pressure_key_ = plist_.get<std::string>("pressure key", "pressure");
  dependencies_.insert(pressure_key_);
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void WRMEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  Epetra_MultiVector& sat_c = *result->ViewComponent("cell", false);
  const Epetra_MultiVector& pres_c = *S->GetFieldData(pressure_key_)->ViewComponent("cell", false);
  const double patm = *S->GetScalarData("atmospheric_pressure");

  int ncells = sat_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    sat_c[0][c] = wrm_->second[(*wrm_->first)[c]]->saturation(patm - pres_c[0][c]);
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void WRMEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  Epetra_MultiVector& sat_c = *result->ViewComponent("cell", false);
  const Epetra_MultiVector& pres_c = *S->GetFieldData(pressure_key_)->ViewComponent("cell", false);
  const double patm = *S->GetScalarData("atmospheric_pressure");

  int ncells = sat_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    // Negative sign indicates that dSdP = -dSdPc.
    sat_c[0][c] = -wrm_->second[(*wrm_->first)[c]]->dSdPc(patm - pres_c[0][c]);
  }
}

}  // namespace Flow
}  // namespace Amanzi
