/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Model for evaluating s_ice, s_liquid, and s_gas.
*/

#ifndef AMANZI_FLOWRELATIONS_WRM_PERMAFROST_MODEL_
#define AMANZI_FLOWRELATIONS_WRM_PERMAFROST_MODEL_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WRMPermafrostModel {

public:
  // required methods from the base class
  virtual void saturations(double pc_liq, double pc_ice, double (&sats)[3]) = 0;
  virtual void dsaturations_dpc_liq(double pc_liq, double pc_ice,
          double (&dsats)[3]) = 0;
  virtual void dsaturations_dpc_ice(double pc_liq, double pc_ice,
          double (&dsats)[3]) = 0;
};

typedef
std::pair<std::string, Teuchos::RCP<WRMPermafrostModel> > WRMPermafrostModelRegionPair;

typedef
std::vector<WRMPermafrostModelRegionPair> WRMPermafrostModelRegionPairList;

} //namespace
} //namespace
} //namespace

#endif
