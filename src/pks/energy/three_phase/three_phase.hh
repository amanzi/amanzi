/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Fully three-phase (air, water, ice) permafrost energy equation.
------------------------------------------------------------------------- */

#ifndef PKS_ENERGY_THREE_PHASE_HH_
#define PKS_ENERGY_THREE_PHASE_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "two_phase.hh"

namespace Amanzi {
namespace Energy {

class ThreePhase : public TwoPhase {

public:
  ThreePhase() {};
  ThreePhase(Teuchos::ParameterList& plist, const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& solution);

protected:
  virtual void SetupPhysicalEvaluators_(const Teuchos::RCP<State>& S);

private:
  static RegisteredPKFactory<ThreePhase> reg_;
  
  friend  class MPCCoupledFlowEnergy;
};

} // namespace Energy
} // namespace Amanzi

#endif
