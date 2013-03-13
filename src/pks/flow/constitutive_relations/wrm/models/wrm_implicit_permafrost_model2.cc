#include <boost/math/tools/roots.hpp>

#include "dbc.hh"
#include "wrm.hh"
#include "wrm_implicit_permafrost_model2.hh"


namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<WRMPermafrostModel,WRMImplicitPermafrostModel2> WRMImplicitPermafrostModel2::factory_("newer permafrost model");

bool WRMImplicitPermafrostModel2::saturations_if_saturated_(double pc_liq,
        double pc_ice, double (&sats)[3]) {
  if (pc_liq <= 0.) { // fully saturated, s_g = 0, s_l/s_i by S_star(pc_ic)
    sats[0] = 0.; // gas
    sats[1] = wrm_->saturation(pc_ice); // liquid
    sats[2] = 1.0 - sats[1];  // ice
    return true;
  }
  return false;
}

bool WRMImplicitPermafrostModel2::saturations_if_above_freezing_(double pc_liq,
        double pc_ice, double (&sats)[3]) {
  if (pc_ice == 0.) {  // above freezing, s_i = 0, s_l/s_g by usual curve
    sats[2] = 0.; // ice
    sats[1] = wrm_->saturation(pc_liq);  // liquid
    sats[0] = 1.0 - sats[1];  // gas
    return true;
  }
  return false;
}

void WRMImplicitPermafrostModel2::saturations(double pc_liq, double pc_ice,
        double (&sats)[3]) {
  if (saturations_if_above_freezing_(pc_liq, pc_ice, sats)) return;
  if (saturations_if_saturated_(pc_liq, pc_ice, sats)) return;

  // solve implicit equation for s_i
  SatIceFunctor_ func(pc_liq, pc_ice, wrm_);
  Tol_ tol(eps_);
  uintmax_t max_it(max_it_);
  double left = 0.;
  double right = 1.;
  std::pair<double,double> result =
      boost::math::tools::toms748_solve(func, left, right, tol, max_it);
  ASSERT(max_it < max_it_);
  sats[2] = result.first;

  // now liq and gas from ice
  sats[1] = (1.0 - sats[2])*wrm_->saturation(pc_liq);
  sats[0] = 1.0 - sats[1] - sats[2];
  return;
}

void WRMImplicitPermafrostModel2::saturations(double pc_liq, double pc_ice,
        double guess, double (&sats)[3]) {
  if (saturations_if_above_freezing_(pc_liq, pc_ice, sats)) return;
  if (saturations_if_saturated_(pc_liq, pc_ice, sats)) return;

  // temperature below 0
  // guess must be away from zero
  if (guess <= 0.) guess = 1.e-8;

  // solve implicit system for sat_i
  SatIceFunctor_ func(pc_liq, pc_ice, wrm_);
  Tol_ tol(eps_);
  double factor = 2.0;
  uintmax_t max_it(max_it_);
  std::pair<double,double> result =
      boost::math::tools::bracket_and_solve_root(func, guess, factor, false, tol, max_it);
  ASSERT(max_it < max_it_);
  sats[2] = result.first;

  // now liq and gas from ice
  sats[1] = (1.0 - sats[2])*wrm_->saturation(pc_liq);
  sats[0] = 1.0 - sats[1] - sats[2];
  return;
};


void WRMImplicitPermafrostModel2::dsaturations_dpc_liq(double pc_liq, double pc_ice,
        double (&dsats)[3]) {
  saturations(pc_liq, pc_ice, dsats);
  dsaturations_dpc_liq(dsats[2], pc_liq, pc_ice, dsats);
};


void WRMImplicitPermafrostModel2::dsaturations_dpc_ice(double pc_liq, double pc_ice,
        double (&dsats)[3]) {
  saturations(pc_liq, pc_ice, dsats);
  dsaturations_dpc_ice(dsats[2], pc_liq, pc_ice, dsats);
};

void WRMImplicitPermafrostModel2::dsaturations_dpc_liq(double s_i, double pc_liq,
        double pc_ice, double (&dsats)[3]) {
  if (pc_ice == 0.) {
    // above freezing
    dsats[2] = 0.;
    dsats[1] = wrm_->d_saturation(pc_liq);
    dsats[0] = -dsats[1];
    return;
  }

  double sstar =  wrm_->saturation(pc_liq);
  double sstarprime = wrm_->d_saturation(pc_liq);
  double tmp = (1.0 - s_i) * sstar;
  double tmpprime = (1.0 - s_i) * sstarprime;

  if (std::abs(1.0 - tmp - s_i) < 1.e-15) {
    // effectively saturated
    dsats[2] = 0.;
    dsats[1] = 0.;
    dsats[0] = 0.;
    return;
  }

  double tmp2 = - wrm_->d_saturation( pc_ice + wrm_->capillaryPressure( tmp + s_i))
      * wrm_->d_capillaryPressure( tmp + s_i );

  double numer = tmpprime * (1 + tmp2);
  double denom = - sstar + tmp2 * (1.0 - sstar);

  dsats[2] = -numer / denom;
  dsats[1] = tmpprime - dsats[2] * sstar;
  dsats[0] = - dsats[1] - dsats[2];
  return;
};


void WRMImplicitPermafrostModel2::dsaturations_dpc_ice(double s_i, double pc_liq,
        double pc_ice, double (&dsats)[3]) {
  if (pc_ice == 0.) {
    // above freezing
    dsats[2] = 0.;
    dsats[1] = 0.;
    dsats[0] = 0.;
    return;
  }

  double sstar =  wrm_->saturation(pc_liq);
  double tmp = (1.0 - s_i) * sstar;

  if (std::abs(1.0 - tmp - s_i) < 1.e-15) {
    // effectively saturated
    dsats[0] = 0.;
    dsats[1] = wrm_->d_saturation(pc_ice);
    dsats[2] = - dsats[1];
    return;
  }

  double tmp2 = - wrm_->d_saturation( pc_ice + wrm_->capillaryPressure( tmp + s_i));

  double numer = tmp2;
  double denom = - sstar + tmp2 * wrm_->d_capillaryPressure(tmp + s_i) * (1.0 - sstar);

  dsats[2] = -numer / denom;
  dsats[1] = - dsats[2]*sstar;
  dsats[0] = - dsats[1] - dsats[2];
};



} // namespace
} // namespace
} // namespace
