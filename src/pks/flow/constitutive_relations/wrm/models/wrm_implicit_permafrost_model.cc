#include <boost/math/tools/roots.hpp>
#include <cmath>

#include "Epetra_SerialDenseMatrix.h"

#include "dbc.hh"
#include "errors.hh"
#include "wrm.hh"
#include "wrm_implicit_permafrost_model.hh"


namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<WRMPermafrostModel,WRMImplicitPermafrostModel> WRMImplicitPermafrostModel::factory_("permafrost model");


// Above freezing calculation methods:
// -- saturation calculation, above freezing
bool WRMImplicitPermafrostModel::sats_unfrozen_(double pc_liq,
        double pc_ice, double (&sats)[3]) {
  if (pc_ice <= 0.) {  // above freezing, s_i = 0, s_l/s_g by usual curve
    sats[2] = 0.; // ice
    sats[1] = wrm_->saturation(pc_liq);  // liquid
    sats[0] = 1.0 - sats[1];  // gas
    return true;
  }
  return false;
}


// -- ds_dpcliq calculation, above freezing
bool WRMImplicitPermafrostModel::dsats_dpc_liq_unfrozen_(double pc_liq,
        double pc_ice, double (&dsats)[3]) {
  if (pc_ice <= 0.) {  // above freezing, s_i = 0, s_l/s_g by usual curve
    dsats[2] = 0.; // ice
    dsats[1] = wrm_->d_saturation(pc_liq);  // liquid
    dsats[0] = - dsats[1];  // gas
    return true;
  }
  return false;
}


// -- ds_dpcice calculation, above freezing
bool WRMImplicitPermafrostModel::dsats_dpc_ice_unfrozen_(double pc_liq,
        double pc_ice, double (&dsats)[3]) {
  if (pc_ice <= 0.) {  // above freezing, s_i = 0, s_l/s_g by usual curve
    dsats[2] = 0.; // ice
    dsats[1] = 0.;
    dsats[0] = 0.;
    return true;
  }
  return false;
}


// saturated calculation methods:
// -- saturation calculation, saturated
bool WRMImplicitPermafrostModel::sats_saturated_(double pc_liq,
        double pc_ice, double (&sats)[3]) {
  if (pc_liq <= 0.) { // fully saturated, s_g = 0, s_l/s_i by S_star(pc_ic)
    sats[0] = 0.; // gas
    sats[1] = wrm_->saturation(pc_ice); // liquid
    sats[2] = 1.0 - sats[1];  // ice
    return true;
  }
  return false;
}


// -- ds_dpcliq calculation, saturated
bool WRMImplicitPermafrostModel::dsats_dpc_liq_saturated_(double pc_liq,
        double pc_ice, double (&dsats)[3]) {
  if (pc_liq <= 0.) { // fully saturated, s_g = 0, s_l/s_i by S_star(pc_ic)
    dsats[0] = 0.; // gas
    dsats[1] = 0.;
    dsats[2] = 0.;
    return true;
  }
  return false;
}


// -- ds_dpcice calculation, saturated
bool WRMImplicitPermafrostModel::dsats_dpc_ice_saturated_(double pc_liq,
        double pc_ice, double (&dsats)[3]) {
  if (pc_liq <= 0.) { // fully saturated, s_g = 0, s_l/s_i by S_star(pc_ic)
    dsats[0] = 0.; // gas
    dsats[1] = wrm_->d_saturation(pc_ice); // liquid
    dsats[2] = - dsats[1];  // ice
    return true;
  }
  return false;
}


// partially frozen, unsaturated calculations
// -- saturation calculation, partially frozen, unsaturated
bool WRMImplicitPermafrostModel::sats_frozen_unsaturated_(double pc_liq,
        double pc_ice, double (&sats)[3]) {
  double si = si_frozen_unsaturated_(pc_liq, pc_ice);
  sats[2] = si;
  sats[1] = (1. - si) * wrm_->saturation(pc_liq);
  sats[0] = 1. - si - sats[1];
  return true;
}

// -- ds_dpcliq calculation, partially frozen, unsaturated
bool WRMImplicitPermafrostModel::dsats_dpc_liq_frozen_unsaturated_(double pc_liq,
        double pc_ice, double (&dsats)[3]) {
  double si = si_frozen_unsaturated_(pc_liq, pc_ice);
  double dsi_dpcliq = dsi_dpc_liq_frozen_unsaturated_(pc_liq, pc_ice, si);
  dsats[2] = dsi_dpcliq;
  dsats[1] = (1. - si) * wrm_->d_saturation(pc_liq) - dsi_dpcliq * wrm_->saturation(pc_liq);
  dsats[0] = - dsats[1] - dsats[2];
  return true;
}

// -- ds_dpcice calculation, partially frozen, unsaturated
bool WRMImplicitPermafrostModel::dsats_dpc_ice_frozen_unsaturated_(double pc_liq,
        double pc_ice, double (&dsats)[3]) {
  double si = si_frozen_unsaturated_(pc_liq, pc_ice);
  double dsi_dpcice = dsi_dpc_ice_frozen_unsaturated_(pc_liq, pc_ice, si);
  dsats[2] = dsi_dpcice;
  dsats[1] = - dsi_dpcice * wrm_->saturation(pc_liq);
  dsats[0] = - dsats[1] - dsats[2];
  return true;
}


// -- si calculation, partially frozen, unsaturated
double WRMImplicitPermafrostModel::si_frozen_unsaturated_(double pc_liq, double pc_ice) {
  // check if we are in the splined region
  double cutoff(0.), si_cutoff(0.);
  DetermineSplineCutoff_(pc_liq, pc_ice, cutoff, si_cutoff);
  if (pc_liq > cutoff) {
    // outside of the spline
    return si_frozen_unsaturated_nospline_(pc_liq, pc_ice);
  } else {
    // fit spline, evaluate
    double spline[3];
    FitSpline_(pc_ice, cutoff, si_cutoff, spline);
    return (spline[0] * pc_liq + spline[1]) * pc_liq + spline[2];
  }
}


// -- dsi_dpcliq calculation, partially frozen, unsaturated
double WRMImplicitPermafrostModel::dsi_dpc_liq_frozen_unsaturated_(double pc_liq,
        double pc_ice, double si) {
  // check if we are in the splined region
  double cutoff(0.), si_cutoff(0.);
  DetermineSplineCutoff_(pc_liq, pc_ice, cutoff, si_cutoff);
  if (pc_liq > cutoff) {
    // outside of the spline
    return dsi_dpc_liq_frozen_unsaturated_nospline_(pc_liq, pc_ice, si);
  } else {
    // fit spline, evaluate
    double spline[3];
    FitSpline_(pc_ice, cutoff, si_cutoff, spline);
    return 2 * spline[0] * pc_liq + spline[1];
  }
}


// -- dsi_dpcice calculation, partially frozen, unsaturated
double WRMImplicitPermafrostModel::dsi_dpc_ice_frozen_unsaturated_(double pc_liq,
        double pc_ice, double si) {
  // check if we are in the splined region
  double cutoff(0.), si_cutoff(0.);
  DetermineSplineCutoff_(pc_liq, pc_ice, cutoff, si_cutoff);
  if (pc_liq > cutoff) {
    // outside of the spline
    return dsi_dpc_ice_frozen_unsaturated_nospline_(pc_liq, pc_ice, si);
  } else {
    // fit two splines, difference neighboring splines
    double spline1[3], spline2[3];
    FitSpline_(pc_ice, cutoff, si_cutoff, spline1);

    double pc_ice2 = pc_ice + 10.;
    double si_cutoff2 = si_frozen_unsaturated_nospline_(cutoff, pc_ice2);
    FitSpline_(pc_ice2, cutoff, si_cutoff2, spline2);

    double dspline[3];
    dspline[0] = (spline2[0] - spline1[0]) / 10.;
    dspline[1] = (spline2[1] - spline1[1]) / 10.;
    dspline[2] = (spline2[2] - spline1[2]) / 10.;

    return (dspline[0] * pc_liq + dspline[1]) * pc_liq + dspline[2];
  }
}


// Helper methods for spline
// -- Determine the point beyond which the spline is not needed
bool WRMImplicitPermafrostModel::DetermineSplineCutoff_(double pc_liq, double pc_ice,
        double& cutoff, double& si) {
  cutoff = std::exp(std::floor(std::log(pc_liq)));
  bool done(false);
  while (!done) {
    si = si_frozen_unsaturated_nospline_(pc_liq, pc_ice);
    double sstar = wrm_->saturation(cutoff);
    if ((1. - si) * sstar + si < (1. - 1.e-12)) {
      done = true;
    } else {
      cutoff = std::exp(std::log(cutoff) + 1.);
    }
  }
  return true;
}


// -- Determine the coefficients of the spline
bool WRMImplicitPermafrostModel::FitSpline_(double pc_ice, double cutoff,
        double si_cutoff, double (&coefs)[3]) {
  double dsi_cutoff = dsi_dpc_liq_frozen_unsaturated_nospline_(cutoff, pc_ice, si_cutoff);


  // spline given by si = a * pc^2 + b * pc + c
  // constraint equations:
  //   1. si(0) = si_saturated  (right-continuous)
  //   2. si(cutoff) = si_cutoff (left-continuous)
  //   3. dsi_dpcliq(cutoff) = dsi_cutoff (continuous left derivative)
  // solve for (a,b,c)

  // c
  sats_saturated_(-1.0, pc_ice, coefs);

  // a,b require the inversion of a linear system of equations
  Epetra_SerialDenseMatrix M(2, 2);
  double rhs[2];

  // -- equation 3
  M(0,0) = 2. * cutoff;
  M(0,1) = 1.;
  rhs[0] = dsi_cutoff;

  // -- equation 2
  M(1,0) = cutoff * cutoff;
  M(1,1) = cutoff;
  rhs[1] = si_cutoff - coefs[2];

  // invert
  double detM = M(0,0)*M(1,1) - M(0,1)*M(1,0);
  coefs[0] = M(1,1) / detM * rhs[0] - M(0,1) / detM * rhs[1];
  coefs[1] = M(0,0) / detM * rhs[1] - M(1,0) / detM * rhs[0];
  return true;
}


// -- si calculation, outside of the splined region
double WRMImplicitPermafrostModel::si_frozen_unsaturated_nospline_(double pc_liq,
        double pc_ice) {
  // solve implicit equation for s_i
  SatIceFunctor_ func(pc_liq, pc_ice, wrm_);
  Tol_ tol(eps_);
  uintmax_t max_it(max_it_);
  double left = 0.;
  double right = 1.;
  std::pair<double,double> result;
  try {
    result =
        boost::math::tools::toms748_solve(func, left, right, tol, max_it);
  } catch (const std::exception& e) {
    std::cout << "WRMImplicitPermafrostModel failed: " << e.what() << std::endl;
    Exceptions::amanzi_throw(Errors::CutTimeStep());
  }

  ASSERT(max_it < max_it_);
  //  std::cout << " took " << max_it << " steps";
  return result.first;
}


// -- dsi_dpcliq calculation, outside of the splined region
double WRMImplicitPermafrostModel::dsi_dpc_liq_frozen_unsaturated_nospline_(double pc_liq,
        double pc_ice, double si) {
  // differentiate the implicit functor, solve for dsi_dpcliq
  double sstar =  wrm_->saturation(pc_liq);
  double sstarprime = wrm_->d_saturation(pc_liq);
  double tmp = (1.0 - si) * sstar;
  double tmpprime = (1.0 - si) * sstarprime;
  double G = - wrm_->d_saturation( pc_ice + wrm_->capillaryPressure( tmp + si))
      * wrm_->d_capillaryPressure( tmp + si );

  double numer = tmpprime * (1 + G);
  double denom = - sstar + G * (1.0 - sstar);

  return -numer / denom;
}


// -- dsi_pcice calculation, outside of the splined region
double WRMImplicitPermafrostModel::dsi_dpc_ice_frozen_unsaturated_nospline_(double pc_liq,
        double pc_ice, double si) {
  // differentiate the implicit functor, solve for dsi_dpcice
  double sstar =  wrm_->saturation(pc_liq);
  double tmp = (1.0 - si) * sstar;
  double G1 = wrm_->d_saturation( pc_ice + wrm_->capillaryPressure( tmp + si));
  double G2 = wrm_->d_capillaryPressure(tmp + si);

  return -G1 / (sstar + G1 * G2 * (1 - sstar));
}


// PUBLIC METHODS
// Calculate the saturation
void WRMImplicitPermafrostModel::saturations(double pc_liq, double pc_ice,
        double (&sats)[3]) {
  if (sats_unfrozen_(pc_liq, pc_ice, sats)) return;
  if (sats_saturated_(pc_liq, pc_ice, sats)) return;
  sats_frozen_unsaturated_(pc_liq, pc_ice, sats);
  return;
}

void WRMImplicitPermafrostModel::dsaturations_dpc_liq(double pc_liq, double pc_ice,
        double (&dsats)[3]) {
  if (dsats_dpc_liq_unfrozen_(pc_liq, pc_ice, dsats)) return;
  if (dsats_dpc_liq_saturated_(pc_liq, pc_ice, dsats)) return;
  dsats_dpc_liq_frozen_unsaturated_(pc_liq, pc_ice, dsats);
};

void WRMImplicitPermafrostModel::dsaturations_dpc_ice(double pc_liq, double pc_ice,
        double (&dsats)[3]) {
  if (dsats_dpc_ice_unfrozen_(pc_liq, pc_ice, dsats)) return;
  if (dsats_dpc_ice_saturated_(pc_liq, pc_ice, dsats)) return;
  dsats_dpc_ice_frozen_unsaturated_(pc_liq, pc_ice, dsats);
};


} // namespace
} // namespace
} // namespace
