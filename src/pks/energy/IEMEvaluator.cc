/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  The internal energy model evaluator simply calls the IEM with 
  the correct arguments.
*/

#include "IEMEvaluator.hh"
#include "IEMFactory.hh"

namespace Amanzi {
namespace Energy {

/* ******************************************************************
* Constructor.
****************************************************************** */
IEMEvaluator::IEMEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  AMANZI_ASSERT(plist_.isSublist("IEM parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("IEM parameters");
  IEMFactory fac;
  iem_ = fac.CreateIEM(sublist);

  InitializeFromPlist_();
}


/* ******************************************************************
* Copy constructor with initialization.
****************************************************************** */
IEMEvaluator::IEMEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<IEM>& iem) :
    SecondaryVariableFieldEvaluator(plist),
    iem_(iem) {
  InitializeFromPlist_();
}


/* ******************************************************************
* Copy constructors.
****************************************************************** */
IEMEvaluator::IEMEvaluator(const IEMEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    iem_(other.iem_),
    temp_key_(other.temp_key_) {};


Teuchos::RCP<FieldEvaluator> IEMEvaluator::Clone() const {
  return Teuchos::rcp(new IEMEvaluator(*this));
}


/* ******************************************************************
* Itialization.
****************************************************************** */
void IEMEvaluator::InitializeFromPlist_()
{
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("internal energy key");
  }

  // Set up my dependencies.
  std::string domain = Keys::getDomain(my_key_);
  temp_key_ = plist_.get<std::string>("temperature key", Keys::getKey(domain, "temperature"));
  dependencies_.insert(temp_key_);
}


/* ******************************************************************
* Evaluation body.
****************************************************************** */
void IEMEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S, const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);

  for (CompositeVector::name_iterator comp = result->begin(); comp != result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp, false);

    int ncomp = result->size(*comp, false);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] = iem_->InternalEnergy(temp_v[0][i]);
    }
  }
}


/* ******************************************************************
* Evaluation of field derivatives.
****************************************************************** */
void IEMEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  AMANZI_ASSERT(wrt_key == temp_key_);
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);

  for (CompositeVector::name_iterator comp = result->begin(); comp != result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp, false);

    int ncomp = result->size(*comp, false);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] = iem_->DInternalEnergyDT(temp_v[0][i]);
    }
  }
}

}  // namespace Energy
}  // namespace Amanzi
