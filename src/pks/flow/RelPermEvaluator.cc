/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Rel perm( pc ( sat ) ).
*/

#include "FlowDefs.hh"
#include "RelPermEvaluator.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Two constructors.
****************************************************************** */
RelPermEvaluator::RelPermEvaluator(Teuchos::ParameterList& plist,
                                   const Teuchos::RCP<WRMPartition>& wrm) :
    SecondaryVariableFieldEvaluator(plist),
    wrm_(wrm),
    min_value_(0.0),
    max_value_(1.0) {
  InitializeFromPlist_();
}

RelPermEvaluator::RelPermEvaluator(const RelPermEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    wrm_(other.wrm_),
    pressure_key_(other.pressure_key_),
    min_value_(other.min_value_),
    max_value_(other.max_value_) {};


/* ******************************************************************
* Copy constructor.
****************************************************************** */
Teuchos::RCP<FieldEvaluator> RelPermEvaluator::Clone() const {
  return Teuchos::rcp(new RelPermEvaluator(*this));
}


/* ******************************************************************
* Initialization.
****************************************************************** */
void RelPermEvaluator::InitializeFromPlist_()
{
  // my keys is for rel perm.
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("rel perm key", "relative_permeability");
  }

  // my dependency is pressure.
  pressure_key_ = plist_.get<std::string>("pressure key", "pressure");
  dependencies_.insert(pressure_key_);
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void RelPermEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result)
{
  // Evaluate the model to calculate krel on cells.
  Epetra_MultiVector& relperm_c = *result->ViewComponent("cell", false);
  const Epetra_MultiVector& pres_c = *S->GetFieldData(pressure_key_)->ViewComponent("cell", false);

  int ncells = relperm_c.MyLength();
  double patm = FLOW_PRESSURE_ATMOSPHERIC;

  for (int c = 0; c != ncells; ++c) {
    relperm_c[0][c] = wrm_->second[(*wrm_->first)[c]]->k_relative(patm - pres_c[0][c]);
  }
}


/* ******************************************************************
* Required member function.
****************************************************************** */
void RelPermEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result)
{
  ASSERT(wrt_key == pressure_key_);

  Epetra_MultiVector& relperm_c = *result->ViewComponent("cell", false);
  const Epetra_MultiVector& pres_c = *S->GetFieldData(pressure_key_)->ViewComponent("cell", false);

  int ncells = relperm_c.MyLength();
  double patm = FLOW_PRESSURE_ATMOSPHERIC;

  for (int c = 0; c != ncells; ++c) {
    relperm_c[0][c] = wrm_->second[(*wrm_->first)[c]]->dKdPc(patm - pres_c[0][c]);
  }
}


/* ******************************************************************
* Single pressure functions.
****************************************************************** */
double RelPermEvaluator::Value(int c, double p) const {
  return wrm_->second[(*wrm_->first)[c]]->k_relative(FLOW_PRESSURE_ATMOSPHERIC - p);
}


double RelPermEvaluator::Derivative(int c, double p) const {
  return wrm_->second[(*wrm_->first)[c]]->dKdPc(FLOW_PRESSURE_ATMOSPHERIC - p);
}


/* ******************************************************************
* CV pressure functions.
****************************************************************** */
void RelPermEvaluator::Value(Teuchos::RCP<CompositeVector>& p, Teuchos::RCP<CompositeVector>& relperm)
{
  Epetra_MultiVector& relperm_c = *relperm->ViewComponent("cell", false);
  const Epetra_MultiVector& p_c = *p->ViewComponent("cell", false);

  int ncells = relperm_c.MyLength();
  double patm = FLOW_PRESSURE_ATMOSPHERIC;

  for (int c = 0; c != ncells; ++c) {
    relperm_c[0][c] = wrm_->second[(*wrm_->first)[c]]->k_relative(patm - p_c[0][c]);
  }
}

}  // namespace Flow
}  // namespace Amanzi
