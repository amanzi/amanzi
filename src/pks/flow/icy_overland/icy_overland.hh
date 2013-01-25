/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#ifndef PK_FLOW_ICY_OVERLAND_HH_
#define PK_FLOW_ICY_OVERLAND_HH_

#include "Teuchos_TimeMonitor.hpp"


#include "pk_factory.hh"
#include "overland.hh"

namespace Amanzi {

namespace Operators { class Upwinding; }

namespace Flow {

class IcyOverlandFlow : public OverlandFlow {

public:
  IcyOverlandFlow(Teuchos::ParameterList& plist,
                  const Teuchos::RCP<TreeVector>& solution) :
      PKDefaultBase(plist, solution),
      OverlandFlow(plist, solution) {}

  // Virtual destructor
  virtual ~IcyOverlandFlow() {}

protected:
  // setup methods
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);

  bool UpdatePermeabilityData_(const Teuchos::Ptr<State>& S);

private:
  // factory registration
  static RegisteredPKFactory<IcyOverlandFlow> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
