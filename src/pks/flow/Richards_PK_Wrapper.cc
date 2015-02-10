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
  glist_ = Teuchos::rcp(new Teuchos::ParameterList(*global_list));
  //We need the flow list
   

  std::string pk_name = pk_tree.name();
  const char* result = pk_name.data();
  while ((result = std::strstr(result, "->")) != NULL) {
    result += 2;
    pk_name = result;   
  }

  // construct
   pk_ = Teuchos::rcp(new Richards_PK(*glist_, pk_name, S_));
}

bool Richards_PK_Wrapper::AdvanceStep(double t_old, double t_new) {
  bool failed = false;
  double dt = t_new - t_old;
  double dt_actual(dt);
  int ierr;
  failed = pk_->Advance(dt, dt_actual);
  if (std::abs(dt - dt_actual) > 1.e-10) {
    failed = true;
  }
  if (failed){
    Teuchos::OSTab tab = pk_->vo_->getOSTab();
    *(pk_->vo_->os()) << "Step failed " << std::endl;
    //std::cout<<"Step failed " << std::endl;
  }

  return failed;
}

}
}
