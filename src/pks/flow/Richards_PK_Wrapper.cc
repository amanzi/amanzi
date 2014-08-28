/*
  License: see $AMANZI_DIR/COPYRIGHT
  Authors: Ethan Coon

  Temporary wrapper converting the Richards_PK, which inherits from 
  BDFFnBase<CompositeVector>, to use TreeVectors.

*/


#include "Richards_PK.hh"
#include "Richards_PK_Wrapper.hh"

namespace Amanzi {
namespace Flow {

Richards_PK_Wrapper::Richards_PK_Wrapper(Teuchos::ParameterList& pk_tree,
        const Teuchos::RCP<Teuchos::ParameterList>& global_list,
        const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& soln) :
    S_(S),
    soln_(soln)
{
  // Richards expects a single global list with sublist Flow
  glist_ = Teuchos::rcp(new ParameterList(*global_list));
  glist_->set("Flow", global_list->sublist("PKs").sublist(pk_tree.name()));

  // construct
  pk_ = Teuchos::rcp(new Richards_PK(*glist_, S_));
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
