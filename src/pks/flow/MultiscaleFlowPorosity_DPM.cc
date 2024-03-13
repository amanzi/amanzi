/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#include <string>

#include "FlowDefs.hh"
#include "MultiscaleFlowPorosity_DPM.hh"
#include "WRMFactory.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* This model is minor extension of the WRM.
****************************************************************** */
MultiscaleFlowPorosity_DPM::MultiscaleFlowPorosity_DPM(Teuchos::ParameterList& plist)
{
  WRMFactory factory;
  wrm_ = factory.Create(plist);

  Teuchos::ParameterList& slist = plist.sublist("dual porosity parameters");
  alpha_ = slist.get<double>("mass transfer coefficient", 0.0);
  tol_ = slist.get<double>("tolerance", FLOW_DPM_NEWTON_TOLERANCE);
  atm_pressure_ = plist.get<double>("atmospheric pressure", FLOW_PRESSURE_ATMOSPHERIC);
}


/* ******************************************************************
* It should be called only once; otherwise, create an evaluator.
****************************************************************** */
double
MultiscaleFlowPorosity_DPM::ComputeField(double phi, double n_l, double prm)
{
  double pc = atm_pressure_ - prm;
  return wrm_->saturation(pc) * phi * n_l;
}


/* ******************************************************************
* Main capability: cell-based Newton solver. It returns water storage,
* pressure in the matrix. max_itrs is input/output parameter.
****************************************************************** */
WhetStone::DenseVector
MultiscaleFlowPorosity_DPM::WaterContentMatrix(double prf0,
                                               WhetStone::DenseVector& prm,
                                               WhetStone::DenseVector& wcm0,
                                               double dt,
                                               double phi,
                                               double n_l,
                                               double mu_l,
                                               int& max_itrs)
{
  double zoom, pmin, pmax, pcm, pcf0;
  pcf0 = atm_pressure_ - prf0;
  pcm = atm_pressure_ - prm(0);

  zoom = std::fabs(pcm) + atm_pressure_;
  pmin = pcm - zoom;
  pmax = pcm + zoom;

  // setup local parameters
  double sat0, alpha_mod;
  sat0 = wcm0(0) / (phi * n_l);
  alpha_mod = alpha_ * dt / (phi * n_l);

  // setup iterative parameters
  double f0, f1, ds, dp, dsdp, guess, result(pcm);
  double delta(1.0e+10), delta1(1.0e+10), delta2(1.0e+10);
  int count(max_itrs);

  while (--count && (fabs(result * tol_) < fabs(delta))) {
    delta2 = delta1;
    delta1 = delta;

    ds = wrm_->saturation(result) - sat0;
    dp = result - pcf0;
    dsdp = wrm_->dSdPc(result);

    f0 = ds - alpha_mod * dp;
    if (f0 == 0.0) break;

    f1 = dsdp - alpha_mod;
    delta = f0 / f1;

    // If the last two steps have not converged, try bisection:
    if (fabs(delta * 2) > fabs(delta2)) {
      delta = (delta > 0) ? (result - pmin) / 2 : (result - pmax) / 2;
    }
    guess = result;
    result -= delta;
    if (result <= pmin) {
      delta = (guess - pmin) / 2;
      result = guess - delta;
      if ((result == pmin) || (result == pmax)) break;

    } else if (result >= pmax) {
      delta = (guess - pmax) / 2;
      result = guess - delta;
      if ((result == pmin) || (result == pmax)) break;
    }

    // update brackets:
    if (delta > 0.0) {
      pmax = guess;
    } else {
      pmin = guess;
    }
  }
  max_itrs -= count - 1;

  WhetStone::DenseVector wcm1(1);
  prm(0) = atm_pressure_ - result;
  wcm1(0) = wrm_->saturation(result) * phi * n_l;

  return wcm1;
}

} // namespace Flow
} // namespace Amanzi
