/*
  License: see $AMANZI_DIR/COPYRIGHT
  Authors: Ethan Coon

  Temporary wrapper converting the Transport_PK, which inherits from 
  BDFFnBase<CompositeVector>, to use TreeVectors.

*/


#include "Transport_PK.hh"
#include "Transport_PK.hh"

#include "Transport_PK_Wrapper.hh"

namespace Amanzi {
namespace Transport {

Transport_PK_Wrapper::Transport_PK_Wrapper(const Teuchos::RCP<Teuchos::ParameterList>& plist,
        Teuchos::ParameterList& glist,
        const Teuchos::RCP<TreeVector>& soln) :
    soln_(soln),
    glist_(glist) {      
  comp_names_ = plist->get<Teuchos::Array<std::string> >("component names").toVector();

  // Transport PK expects a single global Plist with the flow PK as a sublist
  glist_.set("Transport", *plist);
}

void
Transport_PK_Wrapper::SetState(const Teuchos::RCP<State>& S) {
  // can finally construct, as Transport PK expects state in constructor
  pk_ = Teuchos::rcp(new Transport_PK(glist_, S, comp_names_));
  pk_->SetState(S);
}

bool
Transport_PK_Wrapper::Advance(double dt) {
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
