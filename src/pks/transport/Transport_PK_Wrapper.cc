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

Transport_PK_Wrapper::Transport_PK_Wrapper(Teuchos::ParameterList& pk_tree,
        const Teuchos::RCP<Teuchos::ParameterList>& global_list,
        const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& soln) :
    S_(S),
    soln_(soln)
{
  // Transport expects a single global list with sublist Flow
  glist_ = Teuchos::rcp(new Teuchos::ParameterList(*global_list));
  glist_->set("Transport", global_list->sublist("PKs").sublist(pk_tree.name()));

  // grab the component names
  comp_names_ = glist_->sublist("PKs").sublist(pk_tree.name())
      .get<Teuchos::Array<std::string> >("component names").toVector();
  
  // construct
  pk_ = Teuchos::rcp(new Transport_PK(*glist_, S_, comp_names_));
}


/* ******************************************************************
* Wrapper for new MPC policy.
****************************************************************** */
bool Transport_PK_Wrapper::AdvanceStep(double t_old, double t_new) {
  bool failed = false;
  double dt = t_new - t_old;
  double dt_actual(dt);
  int ierr = pk_->Advance(dt, dt_actual);
  if (std::abs(dt - dt_actual) > 1.e-10 || ierr) {
    failed = true;
  }
  return failed;
}

}  // namespace Transport
}  // namespace Amanzi

