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
#include "overland_head.hh"
#include "icy_height_model.hh"

namespace Amanzi {

namespace Operators { class Upwinding; }

namespace Flow {

namespace FlowRelations { class UnfrozenFractionModel; }

class IcyOverlandFlow : public OverlandHeadFlow {

 public:
  IcyOverlandFlow(Teuchos::ParameterList& plist,
                  Teuchos::ParameterList& FElist,
                  const Teuchos::RCP<TreeVector>& solution) :
      PKDefaultBase(plist, FElist, solution),
      OverlandHeadFlow(plist, FElist, solution) {}

  // Virtual destructor
  virtual ~IcyOverlandFlow() {}

 protected:
  // setup methods
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);

  virtual void UpdateBoundaryConditions_(const Teuchos::Ptr<State>& S);

 protected:
  Teuchos::RCP<FlowRelations::IcyHeightModel> icy_height_model_;
  Teuchos::RCP<FlowRelations::UnfrozenFractionModel> uf_model_;

 private:
  // factory registration
  static RegisteredPKFactory<IcyOverlandFlow> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
