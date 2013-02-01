/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  A capillary pressure model based upon something other than p_atm - p.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_PC_LIQ_ATM_
#define AMANZI_FLOW_RELATIONS_PC_LIQ_ATM_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class PCLiqAtm {

public:
  explicit PCLiqAtm(Teuchos::ParameterList& plist) {}

  // required methods from the base class
  virtual double CapillaryPressure(double p, double p_atm) { return p_atm - p; }
  virtual double DCapillaryPressureDp(double p, double p_atm) { return -1.; }
  virtual double DCapillaryPressureDpatm(double p, double p_atm) { return 1.; }
};

} //namespace
} //namespace
} //namespace

#endif
