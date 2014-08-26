/*
  License: see $AMANZI_DIR/COPYRIGHT
  Authors: Ethan Coon

  Temporary wrapper converting the Darcy_PK, which inherits from 
  BDFFnBase<CompositeVector>, to use TreeVectors.

*/


#include "Darcy_PK.hh"
#include "Darcy_PK.hh"

#include "Darcy_PK_Wrapper.hh"

namespace Amanzi {
namespace Flow {

Darcy_PK_Wrapper::Darcy_PK_Wrapper(const Teuchos::RCP<Teuchos::ParameterList>& plist,
        Teuchos::ParameterList& glist,
        const Teuchos::RCP<TreeVector>& soln) :
    soln_(soln),
    glist_(glist) {      
  // Darcy PK expects a single global Plist with the flow PK as a sublist
  glist_.set("Flow", *plist);
}

void
Darcy_PK_Wrapper::SetState(const Teuchos::RCP<State>& S) {
  // can finally construct, as Darcy PK expects state in constructor
  pk_ = Teuchos::rcp(new Darcy_PK(glist_, S));
  pk_->SetState(S);
}

bool
Darcy_PK_Wrapper::Advance(double dt) {
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
