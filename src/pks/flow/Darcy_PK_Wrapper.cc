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

Darcy_PK_Wrapper::Darcy_PK_Wrapper(Teuchos::ParameterList& pk_tree,
        const Teuchos::RCP<Teuchos::ParameterList>& global_list,
        const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& soln) :
    S_(S),
    soln_(soln)
{
  // Darcy expects a single global list with sublist Flow
  glist_ = Teuchos::rcp(new Teuchos::ParameterList(*global_list));
  glist_->set("Flow", global_list->sublist("PKs").sublist(pk_tree.name()));

  // construct
  pk_ = Teuchos::rcp(new Darcy_PK(*glist_, S_));
}


bool
Darcy_PK_Wrapper::AdvanceStep(double t_old, double t_new) {
  bool failed = false;
  double dt = t_new - t_old;
  double dt_actual(dt);
  int ierr = pk_->Advance(dt, dt_actual);
  if (std::abs(dt - dt_actual) > 1.e-10 || ierr) {
    failed = true;
  }
  return failed;
}

} // namespace Flow
} // namespace Amanzi

