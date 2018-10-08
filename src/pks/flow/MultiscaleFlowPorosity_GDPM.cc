/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>

#include "FlowDefs.hh"
#include "MultiscaleFlowPorosity_GDPM.hh"
#include "WRMFactory.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* This model is minor extension of the WRM.
****************************************************************** */
MultiscaleFlowPorosity_GDPM::MultiscaleFlowPorosity_GDPM(Teuchos::ParameterList& plist)
{
  WRMFactory factory;
  wrm_ = factory.Create(plist);

  auto& sublist = plist.sublist("generalized dual porosity parameters");
  matrix_nodes_ = sublist.get<int>("number of matrix nodes");

  // depth is defined for each matrix block as A_m / V_m, so in general,
  // it depends on geometry 
  depth_ = sublist.get<double>("matrix depth");
  tau_ = sublist.get<double>("matrix tortuosity");
  double Ka = sublist.get<double>("absolute permeability");

  tol_ = plist.get<double>("tolerance", FLOW_DPM_NEWTON_TOLERANCE);

  // make uniform mesh inside the matrix
  auto mesh = std::make_shared<WhetStone::DenseVector>(WhetStone::DenseVector(matrix_nodes_ + 1));
  double h = depth_ / matrix_nodes_;
  for (int i = 0; i < matrix_nodes_ + 1; ++i) (*mesh)(i) = h * i;

  // initialize diffusion operator
  auto kr_ = std::make_shared<WhetStone::DenseVector>(WhetStone::DenseVector(matrix_nodes_));
  auto dkdp_ = std::make_shared<WhetStone::DenseVector>(WhetStone::DenseVector(matrix_nodes_));

  op_diff_.Init(mesh);
  op_diff_.Setup(Ka);
  op_diff_.Setup(kr_, dkdp_);
}


/* ******************************************************************
* Compute water content.
****************************************************************** */
double MultiscaleFlowPorosity_GDPM::ComputeField(double phi, double n_l, double pcm)
{
  return wrm_->saturation(pcm) * phi * n_l;
}


/* ******************************************************************
* Main capability: cell-based Newton solver. It returns water content, 
* pressure in the matrix. max_itrs is input/output parameter.
****************************************************************** */
double MultiscaleFlowPorosity_GDPM::WaterContentMatrix(
    double dt, double phi, double n_l, double wcm0, double pcf0, double& pcm, int& max_itrs)
{
  double patm(1e+5), zoom, pmin, pmax;
  zoom = fabs(pcm) + patm;
  pmin = pcm - zoom; 
  pmax = pcm + zoom; 

  // setup local parameters 
  double sat0, alpha_mod, alpha_(1.0);
  sat0 = wcm0 / (phi * n_l);
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

  pcm = result;
  return wrm_->saturation(pcm) * phi * n_l;
}

}  // namespace Flow
}  // namespace Amanzi
 
