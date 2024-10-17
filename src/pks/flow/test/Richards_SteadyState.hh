/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  This is the flow component of the Amanzi code.

*/

#include <algorithm>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "Richards_PK.hh"

namespace Amanzi {
namespace Flow {

class TI_Specs {
 public:
  double T0, dT0;
  double T1;
  int num_itrs, max_itrs;
};


int
AdvanceToSteadyState(Teuchos::RCP<State> S,
                     Richards_PK& RPK,
                     TI_Specs& ti_specs,
                     Teuchos::RCP<TreeVector> soln)
{
  bool last_step = false;

  std::string passwd("");

  int max_itrs = ti_specs.max_itrs;
  double T0 = ti_specs.T0;
  double T1 = ti_specs.T1;
  double dT0 = ti_specs.dT0;

  double T_physics = 0.0;
  double dT(dT0), dTnext;

  int itrs = 0;
  while (itrs < max_itrs && T_physics < T1) {
    if (itrs == 0) { // initialization of BDF1
      Teuchos::RCP<TreeVector> udot = Teuchos::rcp(new TreeVector(*soln));
      udot->PutScalar(0.0);
      RPK.get_bdf1_dae()->SetInitialState(T0, soln, udot);

      RPK.UpdatePreconditioner(T0, soln, dT0);
    }

    while (RPK.get_bdf1_dae()->AdvanceStep(dT, dTnext, soln)) { dT = dTnext; }
    RPK.VV_ReportSeepageOutflow(S.ptr(), dT);
    RPK.get_bdf1_dae()->CommitSolution(dT, soln);

    T_physics = RPK.get_bdf1_dae()->time();
    dT = dTnext;
    itrs++;

    double Tdiff = T1 - T_physics;
    if (dTnext > Tdiff) {
      dT = Tdiff * 0.99999991; // To avoid hitting the wrong BC
      last_step = true;
    }
    if (last_step && dT < 1e-3) break;

    // reset primary field
    auto pressure_eval =
      Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CompositeVector, CompositeVectorSpace>>(
        S->GetEvaluatorPtr("pressure", Tags::DEFAULT));
    S->GetW<CompositeVector>("pressure", passwd) = *soln->Data();
    pressure_eval->SetChanged();

    // update and swap saturations
    S->GetEvaluator("saturation_liquid").Update(*S, "flow");
    const auto& s_l = S->Get<CompositeVector>("saturation_liquid");
    auto& s_l_prev = S->GetW<CompositeVector>("prev_saturation_liquid", passwd);
    s_l_prev = s_l;

    // update and swap water storage
    S->GetEvaluator("water_storage").Update(*S, "flow");
    const auto& wc = S->Get<CompositeVector>("water_storage");
    auto& wc_prev = S->GetW<CompositeVector>("prev_water_storage", passwd);
    wc_prev = wc;

    // update and swap multiscale water storage
    if (S->HasRecord("water_storage_msp")) {
      const auto& wcm = S->Get<CompositeVector>("water_storage_msp");
      auto& wcm_prev = S->GetW<CompositeVector>("prev_water_storage_msp", passwd);
      wcm_prev = wcm;
    }

    // commit step
    RPK.CommitStep(T_physics - dT, T_physics, Tags::DEFAULT);
    RPK.VV_ReportMultiscale();
  }

  ti_specs.num_itrs = itrs;
  return 0;
}

} // namespace Flow
} // namespace Amanzi
