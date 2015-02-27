/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)

  The WRM Evaluator simply calls the WRM with the correct arguments.
*/

#include "WRM_evaluator.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Two constructors.
****************************************************************** */
WRMEvaluator::WRMEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariablesFieldEvaluator(plist)
{
  ASSERT(plist_.isSublist("WRM parameters"));
  Teuchos::ParameterList wrm_plist = plist_.sublist("WRM parameters");
  // wrm_ = createWRMPartition(wrm_plist);

  InitializeFromPlist_();
}


WRMEvaluator::WRMEvaluator(const WRMEvaluator& other) :
    SecondaryVariablesFieldEvaluator(other),
    capillary_pressure_key_(other.capillary_pressure_key_),
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
  // my key is for saturation
  my_keys_.push_back(plist_.get<std::string>("saturation key", "saturation_liquid"));

  // my dependency is capillary pressure.
  capillary_pressure_key_ = plist_.get<std::string>("capillary pressure key", "capillary_pressure_gas_liquid");
  dependencies_.insert(capillary_pressure_key_);
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void WRMEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const std::vector<Teuchos::Ptr<CompositeVector> >& results)
{
  Epetra_MultiVector& sat_c = *results[0]->ViewComponent("cell", false);
  const Epetra_MultiVector& pres_c = *S->GetFieldData(capillary_pressure_key_)->ViewComponent("cell", false);

  int ncells = sat_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    sat_c[0][c] = wrm_[(*cell2wrm_)[c]]->saturation(pres_c[0][c]);
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void WRMEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key,
    const std::vector<Teuchos::Ptr<CompositeVector> > & results)
{
  Epetra_MultiVector& sat_c = *results[0]->ViewComponent("cell", false);
  const Epetra_MultiVector& pres_c = *S->GetFieldData(capillary_pressure_key_)->ViewComponent("cell", false);

  int ncells = sat_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    sat_c[0][c] = wrm_[(*cell2wrm_)[c]]->dSdPc(pres_c[0][c]);
  }
}

}  // namespace Flow
}  // namespace Amanzi
