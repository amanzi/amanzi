/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  The WRM Evaluator simply calls the WRM with the correct arguments.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "wrm_evaluator.hh"
#include "wrm_factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

WRMEvaluator::WRMEvaluator(Teuchos::ParameterList& wrm_plist) :
    SecondaryVariablesFieldEvaluator(),
    wrm_plist_(wrm_plist) {
  ASSERT(wrm_plist.isSublist("WRM parameters"));
  Teuchos::ParameterList sublist = wrm_plist.sublist("WRM parameters");
  WRMFactory fac;
  wrm_ = fac.createWRM(sublist);
}

WRMEvaluator::WRMEvaluator(Teuchos::ParameterList& wrm_plist, const Teuchos::RCP<WRM>& wrm) :
    SecondaryVariablesFieldEvaluator(),
    wrm_plist_(wrm_plist),
    wrm_(wrm) {}

WRMEvaluator::WRMEvaluator(const WRMEvaluator& other) :
    SecondaryVariablesFieldEvaluator(other),
    wrm_plist_(other.wrm_plist_),
    wrm_(other.wrm_) {}

} //namespace
} //namespace
} //namespace
