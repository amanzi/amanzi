/*
  License: see $AMANZI_DIR/COPYRIGHT
  Authors: Daniil Svyatskiy

  Dummy PK which demonstrates the require interface for PK
  BDFFnBase<CompositeVector>, to use TreeVectors.
*/


#include "Dummy_PK.hh"
#include <iostream>

namespace Amanzi {

Dummy_PK::Dummy_PK(Teuchos::ParameterList& pk_tree,
        const Teuchos::RCP<Teuchos::ParameterList>& global_list,
        const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& soln) :
    soln_(soln)
{
  S_ = S;
  glist_ = Teuchos::rcp(new Teuchos::ParameterList(*global_list));
  ti_list_ = glist_->sublist("cycle driver").sublist("time intervals").sublist("TI 0");
}


bool Dummy_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool failed = false;

  if ((step_count + 2)%3 == 0) {
    failed = true;
    dummy_dt = 0.8*dummy_dt;
    std::cout<<"Step failed\n";
  }
  else {
    failed = false;
    dummy_dt = 1.2*dummy_dt;
    std::cout<<"Step succeed. New time "<<t_new<<"\n";
  }
  
  step_count++;

  return failed;
}

}
