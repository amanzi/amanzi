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
    wrm_plist_(wrm_plist) {
  WRMFactory fac;
  wrm_ = fac.createWRM(wrm_plist_);
}

WRMModel::WRMModel(const WRMModel& other) :
    SecondaryVariablesFieldModel(other),
    wrm_plist_(other.wrm_plist_),
    wrm_(other.wrm_) {}

} //namespace
} //namespace
} //namespace
