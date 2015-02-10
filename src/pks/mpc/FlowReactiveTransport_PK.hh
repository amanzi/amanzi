/*
  License: see $AMANZI_DIR/COPYRIGHT
  Authors: Daniil Svyatskiy

  PK for coupling of Flow PK with Transport_PK and Chemestry_PK

*/


#ifndef AMANZI_FLOWREACTIVETRANSPORT_PK_HH_
#define AMANZI_FLOWREACTIVETRANSPORT_PK_HH_

#include "Teuchos_RCP.hpp"

#include "PK.hh"
#include "PK_Factory.hh"
#include "MPCSubcycled.hh"

namespace Amanzi {

class FlowReactiveTransport_PK : public MPCSubcycled{

public:
  FlowReactiveTransport_PK(Teuchos::ParameterList& pk_tree,
               const Teuchos::RCP<Teuchos::ParameterList>& global_list,
               const Teuchos::RCP<State>& S,
               const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();
  virtual void set_dt(double dt);

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new);

  virtual void CommitStep(double t_old, double t_new);

  std::string name() { return "flow reactive transport";} 

private:

  // factory registration
  static RegisteredPKFactory<FlowReactiveTransport_PK> reg_;
};

} // close namespace Amanzi
#endif
