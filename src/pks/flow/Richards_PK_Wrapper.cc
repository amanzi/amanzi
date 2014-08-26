/*
  License: see $AMANZI_DIR/COPYRIGHT
  Authors: Ethan Coon

  Temporary wrapper converting the Richards_PK, which inherits from 
  BDFFnBase<CompositeVector>, to use TreeVectors.

*/


#include "Richards_PK.hh"
#include "Darcy_PK.hh"

#include "Richards_PK_Wrapper.hh"

namespace Amanzi {
namespace Flow {

Richards_PK_Wrapper::Richards_PK_Wrapper(const Teuchos::RCP<Teuchos::ParameterList>& plist,
        Teuchos::ParameterList& glist,
        const Teuchos::RCP<TreeVector>& soln) :
    soln_(soln),
    glist_(glist) {      
  // Richards PK expects a single global Plist with the flow PK as a sublist
  glist_.set("Flow", *plist);
}

void
Richards_PK_Wrapper::SetState(const Teuchos::RCP<State>& S) {
  // can finally construct, as Richards PK expects state in constructor
  pk_ = Teuchos::rcp(new Richards_PK(glist_, S));
  pk_->SetState(S);
}

bool
Richards_PK_Wrapper::Advance(double dt) {
  bool failed = false;
  double dt_actual(dt);
  int ierr = pk_->Advance(dt, dt_actual);
  if (std::abs(dt - dt_actual) > 1.e-10 || ierr) {
    failed = true;
  }
  return failed;
}



}
}
