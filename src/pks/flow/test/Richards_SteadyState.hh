/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "LinearOperatorFactory.hh"
#include "Richards_PK.hh"

namespace Amanzi {
namespace Flow {

class TI_Specs {
 public:
  double T0, dT0;
  double T1;
  int num_itrs, max_itrs;
};


int AdvanceToSteadyState(
    Richards_PK& RPK, TI_Specs& ti_specs, Teuchos::RCP<CompositeVector> solution)
{
  bool last_step = false;

  int max_itrs = ti_specs.max_itrs;
  double T0 = ti_specs.T0;
  double T1 = ti_specs.T1;
  double dT0 = ti_specs.dT0;

  double T_physics = 0.0;
  double dT = dT0, dTnext;

  int itrs = 0;
  while (itrs < max_itrs && T_physics < T1) {
    if (itrs == 0) {  // initialization of BDF1
      Teuchos::RCP<CompositeVector> udot = Teuchos::rcp(new CompositeVector(*solution));
      udot->PutScalar(0.0);
      RPK.get_bdf1_dae()->SetInitialState(T0, solution, udot);

      RPK.UpdatePreconditioner(T0, solution, dT0);
    }

    while (RPK.get_bdf1_dae()->TimeStep(dT, dTnext, solution)) {
      dT = dTnext;
    }
    RPK.get_bdf1_dae()->CommitSolution(dT, solution);

    T_physics = RPK.get_bdf1_dae()->time();
    dT = dTnext;
    itrs++;

    double Tdiff = T1 - T_physics;
    if (dTnext > Tdiff) {
      dT = Tdiff * 0.99999991;  // To avoid hitting the wrong BC
      last_step = true;
    }
    if (last_step && dT < 1e-3) break;
  }

  ti_specs.num_itrs = itrs;
  return 0;
}

}  // namespace Flow
}  // namespace Amanzi

