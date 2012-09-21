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
  ASSERT(wrm_plist_.isSublist("WRM parameters"));

  WRMFactory fac;

  Teuchos::ParameterList region_list = wrm_plist_.sublist("WRM parameters");
  wrms_ = Teuchos::rcp(new WRMRegionPairList());

  for (Teuchos::ParameterList::ConstIterator lcv=region_list.begin();
       lcv!=region_list.end(); ++lcv) {
    std::string name = lcv->first;
    if (region_list.isSublist(name)) {
      Teuchos::ParameterList sublist = region_list.sublist(name);
      std::string region = sublist.get<std::string>("region");
      wrms_->push_back(std::make_pair(region, fac.createWRM(sublist)));
    } else {
      ASSERT(0);
    }
  }
}

WRMEvaluator::WRMEvaluator(Teuchos::ParameterList& wrm_plist,
                           const Teuchos::RCP<std::vector<WRMRegionPair> >& wrms) :
    SecondaryVariablesFieldEvaluator(),
    wrm_plist_(wrm_plist),
    wrms_(wrms) {}

WRMEvaluator::WRMEvaluator(const WRMEvaluator& other) :
    SecondaryVariablesFieldEvaluator(other),
    wrm_plist_(other.wrm_plist_),
    wrms_(other.wrms_) {}

} //namespace
} //namespace
} //namespace
