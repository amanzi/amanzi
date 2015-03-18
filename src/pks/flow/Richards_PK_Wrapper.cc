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

  // swap saturations (may go to a high-level PK)
  S_->GetFieldEvaluator("saturation_liquid")->HasFieldChanged(S_.ptr(), "flow");
  const CompositeVector& sat = *S_->GetFieldData("saturation_liquid");
  CompositeVector& sat_prev = *S_->GetFieldData("prev_saturation_liquid", passwd);

  CompositeVector sat_prev_copy(sat_prev);
  sat_prev = sat;

  // swap water_content (may go to a high-level PK)
  S_->GetFieldEvaluator("water_content")->HasFieldChanged(S_.ptr(), "flow");
  CompositeVector& wc = *S_->GetFieldData("water_content", "water_content");
  CompositeVector& wc_prev = *S_->GetFieldData("prev_water_content", passwd);

  CompositeVector wc_prev_copy(wc_prev);
  wc_prev = wc;

  bool failed = false;
  failed = pk_->Advance(dt, dt_actual);

  if (std::abs(dt - dt_actual) > 1.e-10) {
    failed = true;
  }
  if (failed) {
    // revover the original primary solution, pressure
    *S_->GetFieldData("pressure", passwd) = pressure_copy;
    pk_->pressure_eval()->SetFieldAsChanged(S_.ptr());

    // revover the original fields
    *S_->GetFieldData("prev_saturation_liquid", passwd) = sat_prev_copy;
    *S_->GetFieldData("prev_water_content", passwd) = wc_prev_copy;

    Teuchos::OSTab tab = pk_->vo_->getOSTab();
    *(pk_->vo_->os()) << "Step failed. Restored pressure, prev_saturation_liquid, prev_water_content" << std::endl;
  }

  return failed;
}

}  // namespace Flow
}  // namespace Amanzi
