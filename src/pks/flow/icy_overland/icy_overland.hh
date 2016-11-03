/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#ifndef PK_FLOW_ICY_OVERLAND_HH_
#define PK_FLOW_ICY_OVERLAND_HH_

#include "Teuchos_TimeMonitor.hpp"

#include "PK_Factory.hh"
#include "overland_pressure.hh"
#include "icy_height_model.hh"

namespace Amanzi {

namespace Operators { class Upwinding; }

namespace Flow {

namespace FlowRelations { class UnfrozenFractionModel; }

class IcyOverlandFlow : public OverlandPressureFlow {

 public:
  IcyOverlandFlow(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& solution) :
    //PKDefaultBase(global_list, pk_tree, solution),
    //PK_Default(global_list, pk_tree, solution),
    PK(pk_tree, global_list, S, solution),
    PK_BDF_Default(pk_tree, global_list, S, solution),
    OverlandPressureFlow(pk_tree, global_list, S, solution) {}

  // Virtual destructor
  virtual ~IcyOverlandFlow() {}

 protected:
  // setup methods
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);

 private:
  // factory registration
  static RegisteredPKFactory<IcyOverlandFlow> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
