/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  The WRM Evaluator simply calls the WRM with the correct arguments.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "iem_evaluator.hh"
#include "iem_factory.hh"

namespace Amanzi {
namespace Energy {
namespace EnergyRelations {

IEMEvaluator::IEMEvaluator(Teuchos::ParameterList& iem_plist) :
    SecondaryVariableFieldEvaluator(),
    iem_plist_(iem_plist) {
  ASSERT(iem_plist.isSublist("IEM parameters"));
  Teuchos::ParameterList sublist = iem_plist.sublist("IEM parameters");
  IEMFactory fac;
  iem_ = fac.createIEM(sublist);

  InitializeFromPlist_();
}

IEMEvaluator::IEMEvaluator(Teuchos::ParameterList& iem_plist, const Teuchos::RCP<IEM>& iem) :
    SecondaryVariableFieldEvaluator(),
    iem_plist_(iem_plist),
    iem_(iem) {
  InitializeFromPlist_();
}

IEMEvaluator::IEMEvaluator(const IEMEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    iem_plist_(other.iem_plist_),
    iem_(other.iem_),
    temp_key_(other.temp_key_) {}

Teuchos::RCP<FieldEvaluator>
IEMEvaluator::Clone() const {
  return Teuchos::rcp(new IEMEvaluator(*this));
}



} //namespace
} //namespace
} //namespace
