/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Fully three-phase (air, water, ice) permafrost energy equation, with only
mobile water.

Inherits TwoPhase instead of EnergyBase to pick up the enthalpy from TwoPhase.
------------------------------------------------------------------------- */

#ifndef PKS_ENERGY_THREE_PHASE_HH_
#define PKS_ENERGY_THREE_PHASE_HH_

#include "pk_factory.hh"
#include "three_phase.hh"

namespace Amanzi {
namespace Energy {

class Interfrost : public ThreePhase {

public:
<<<<<<< HEAD
  Interfrost(Teuchos::Ptr<State> S, const Teuchos::RCP<Teuchos::ParameterList>& plist,
             Teuchos::ParameterList& FElist,
             const Teuchos::RCP<TreeVector>& solution) :
    PKDefaultBase(S, plist, FElist, solution),
    ThreePhase(S, plist, FElist, solution) {}
=======
  Interfrost(const Teuchos::RCP<Teuchos::ParameterList>& plist,
             Teuchos::ParameterList& FElist,
             const Teuchos::RCP<TreeVector>& solution) :
      PKDefaultBase(plist, FElist, solution),
      ThreePhase(plist, FElist, solution) {}
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e

  // Virtual destructor
  virtual ~Interfrost() {}
  
 protected:
  virtual void AddAccumulation_(const Teuchos::Ptr<CompositeVector>& f);
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);

private:
  static RegisteredPKFactory<Interfrost> reg_;

  friend  class MPCCoupledFlowEnergy;
};

} // namespace Energy
} // namespace Amanzi

#endif
