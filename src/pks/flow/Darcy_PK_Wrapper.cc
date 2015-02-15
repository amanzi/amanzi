/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

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
  std::string pk_name = pk_tree.name();
  const char* result = pk_name.data();

  while ((result = std::strstr(result, "->")) != NULL) {
    result += 2;
    pk_name = result;
    
  }
  
  // Darcy expects a single global list with sublist Flow
  pk_ = Teuchos::rcp(new Darcy_PK(global_list, pk_name,  S_));
}


bool
Darcy_PK_Wrapper::AdvanceStep(double t_old, double t_new) {
  bool failed = false;
  double dt = t_new - t_old;
  double dt_actual(dt);
  int ierr;
  failed = pk_->Advance(dt, dt_actual);
  if (std::abs(dt - dt_actual) > 1.e-10) {
    failed = true;
  }
  return failed;
}

}  // namespace Flow
}  // namespace Amanzi

