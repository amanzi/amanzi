/*
  License: see $AMANZI_DIR/COPYRIGHT
  Authors: Ethan Coon

  Temporary wrapper converting the Transport_PK, which inherits from 
  BDFFnBase<CompositeVector>, to use TreeVectors.

*/


#include "Transport_PK.hh"
#include "Transport_PK_Wrapper.hh"

namespace Amanzi {
namespace Transport {

Transport_PK_Wrapper::Transport_PK_Wrapper(
        Teuchos::ParameterList& pk_tree,
        const Teuchos::RCP<Teuchos::ParameterList>& glist,
        const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& soln) :
    glist_(glist),
    S_(S),
    soln_(soln)
{
  std::string pk_name = pk_tree.name();
  const char* result = pk_name.data();
  while ((result = std::strstr(result, "->")) != NULL) {
    result += 2;
    pk_name = result;   
  }

  if (glist_->isSublist("Cycle Driver")){
    if (glist_->sublist("Cycle Driver").isParameter("component names")){
      // grab the component names
      comp_names_ = glist_->sublist("Cycle Driver").get<Teuchos::Array<std::string> >("component names").toVector();
    }else{
      Errors::Message msg("Transport PK: parameter component names is missing.");
      Exceptions::amanzi_throw(msg);
    }
  }else{
    Errors::Message msg("Transport PK: sublist Cycle Driver is missing.");
      Exceptions::amanzi_throw(msg);
  }


  Teuchos::RCP<Teuchos::ParameterList> pk_list =  Teuchos::sublist(glist_, "PKs", true);
  Teuchos::RCP<Teuchos::ParameterList> pk_transp_list =  Teuchos::sublist(pk_list, pk_name, true);

  transport_subcycling = pk_transp_list->get<bool>("transport subcycling", true);

   
  // construct
  pk_ = Teuchos::rcp(new Transport_PK(glist_, S_, pk_name, comp_names_));
}


bool
Transport_PK_Wrapper::AdvanceStep(double t_old, double t_new) {

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

