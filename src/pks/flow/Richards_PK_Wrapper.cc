/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon

  Temporary wrapper converting the Richards_PK, which inherits from 
  BDFFnBase<CompositeVector>, to use TreeVectors.
*/


#include "Richards_PK.hh"
#include "Richards_PK_Wrapper.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Standard PK constructor with four arguments.
****************************************************************** */
Richards_PK_Wrapper::Richards_PK_Wrapper(Teuchos::ParameterList& pk_tree,
        const Teuchos::RCP<Teuchos::ParameterList>& global_list,
        const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& soln) :
    S_(S),
    soln_(soln)
{
  std::string pk_name = pk_tree.name();
  const char* result = pk_name.data();
  while ((result = std::strstr(result, "->")) != NULL) {
    result += 2;
    pk_name = result;   
  }

  pk_ = Teuchos::rcp(new Richards_PK(global_list, pk_name, S_, soln_));
}


/* ******************************************************************
* Advance step with save control of primary field pressure.
****************************************************************** */
bool Richards_PK_Wrapper::AdvanceStep(double t_old, double t_new)
{
  double dt = t_new - t_old;
  double dt_actual(dt);

  // save a copy of pressure
  std::string passwd = pk_->name();
  CompositeVector pressure_copy(*S_->GetFieldData("pressure", passwd));

  bool failed = false;
  failed = pk_->Advance(dt, dt_actual);

  if (std::abs(dt - dt_actual) > 1.e-10) {
    failed = true;
  }
  if (failed) {
    // revover the original solution
    *S_->GetFieldData("pressure", passwd) = pressure_copy;
    pk_->pressure_eval()->SetFieldAsChanged(S_.ptr());

    Teuchos::OSTab tab = pk_->vo_->getOSTab();
    *(pk_->vo_->os()) << "Step failed " << std::endl;
  }

  return failed;
}

}  // namespace Flow
}  // namespace Amanzi
