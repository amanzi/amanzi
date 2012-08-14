/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  The WRM Model simply calls the WRM with the correct arguments.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "wrm_model.hh"
#include "wrm_factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

WRMModel::WRMModel(Teuchos::ParameterList& wrm_plist) :
    SecondaryVariablesFieldModel(),
    wrm_plist_(wrm_plist) {
  ASSERT(wrm_plist.isSublist("WRM parameters"));
  Teuchos::ParameterList sublist = wrm_plist.sublist("WRM parameters");
  WRMFactory fac;
  wrm_ = fac.createWRM(sublist);
}

WRMModel::WRMModel(Teuchos::ParameterList& wrm_plist, const Teuchos::RCP<WRM>& wrm) :
    SecondaryVariablesFieldModel(),
    wrm_plist_(wrm_plist),
    wrm_(wrm) {}

WRMModel::WRMModel(const WRMModel& other) :
    SecondaryVariablesFieldModel(other),
    wrm_plist_(other.wrm_plist_),
    wrm_(other.wrm_) {}

} //namespace
} //namespace
} //namespace
