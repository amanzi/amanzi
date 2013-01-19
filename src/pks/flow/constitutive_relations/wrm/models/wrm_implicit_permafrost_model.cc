#include <boost/math/tools/roots.hpp>

#include "dbc.hh"
#include "wrm.hh"
#include "wrm_implicit_permafrost_model.hh"


namespace Amanzi {
namespace Flow {
namespace FlowRelations {


bool WRMImplicitPermafrostModel::saturations_if_above_freezing_(double pc_liq,
        double pc_ice, double (&sats)[3]) {
  if (pc_ice == 0.) {
    sats[2] = 0.;
    sats[1] = wrm_->saturation(pc_liq);
    sats[0] = 1.0 - sats[1];
    return true;
  } else {
    return false;
  }
}

bool WRMImplicitPermafrostModel::dsaturations_dpc_ice_if_above_freezing_(double pc_liq,
        double pc_ice, double (&dsats)[3]) {
  if (pc_ice == 0.) {
    dsats[2] = 0.;
    dsats[1] = 0.;
    dsats[0] = 0.;
    return true;
  } else {
    return false;
  }
}


bool WRMImplicitPermafrostModel::dsaturations_dpc_liq_if_above_freezing_(double pc_liq,
        double pc_ice, double (&dsats)[3]) {
  if (pc_ice == 0.) {
    dsats[2] = 0.;
    dsats[1] = wrm_->d_saturation(pc_liq);
    dsats[0] = -dsats[1];
    return true;
  } else {
    return false;
  }
}



void WRMImplicitPermafrostModel::saturations(double pc_liq, double pc_ice,
        double (&sats)[3]) {
  return saturations(pc_liq, pc_ice, 0.5, sats);
}

void WRMImplicitPermafrostModel::saturations(double pc_liq, double pc_ice,
        double guess, double (&sats)[3]) {

  if (!saturations_if_above_freezing_(pc_liq, pc_ice, sats)) {
    // temperature below 0
    SatIceFunctor_ func(pc_liq, pc_ice, wrm_);
    Tol_ tol(eps_);

    if (guess <= 0.) guess = 1.e-5;

    double factor = 2.0;
    uintmax_t max_it(max_it_);
    std::pair<double,double> result =
        boost::math::tools::bracket_and_solve_root(func, guess, factor, true, tol, max_it);

    ASSERT(max_it < max_it_);
    sats[2] = result.first;

    // now liq and gas from ice
    sats[1] = (1.0 - sats[2])*wrm_->saturation(pc_liq);
    sats[0] = 1.0 - sats[1] - sats[2];
  }
};


void WRMImplicitPermafrostModel::dsaturations_dpc_liq(double pc_liq, double pc_ice,
        double (&dsats)[3]) {
  if (!dsaturations_dpc_liq_if_above_freezing_(pc_liq, pc_ice, dsats)) {
    saturations(pc_liq, pc_ice, dsats);
    double s_i = dsats[2];
    dsaturations_dpc_liq(s_i, pc_liq, pc_ice, dsats);
  }
};


void WRMImplicitPermafrostModel::dsaturations_dpc_ice(double pc_liq, double pc_ice,
        double (&dsats)[3]) {
  if (!dsaturations_dpc_ice_if_above_freezing_(pc_liq, pc_ice, dsats)) {
    saturations(pc_liq, pc_ice, dsats);
    dsaturations_dpc_ice(dsats[2], pc_liq, pc_ice, dsats);
  }
};

void WRMImplicitPermafrostModel::dsaturations_dpc_liq(double s_i, double pc_liq,
        double pc_ice, double (&dsats)[3]) {
  if (!dsaturations_dpc_liq_if_above_freezing_(pc_liq, pc_ice, dsats)) {
    // finite difference to get a guess
    saturations(pc_liq + 1.e-3, pc_ice, s_i, dsats);
    double guess = (dsats[2] - s_i) / 1.e-3;
    dsaturations_dpc_liq(s_i, pc_liq, pc_ice, guess, dsats);
  }
};

void WRMImplicitPermafrostModel::dsaturations_dpc_ice(double s_i, double pc_liq,
        double pc_ice, double (&dsats)[3]) {
  if (!dsaturations_dpc_ice_if_above_freezing_(pc_liq, pc_ice, dsats)) {
    // finite difference to get a guess
    saturations(pc_liq, pc_ice + 1.e-3, s_i, dsats);
    double guess = (dsats[2] - s_i) / 1.e-3;
    dsaturations_dpc_ice(s_i, pc_liq, pc_ice, guess, dsats);
  }
};

void WRMImplicitPermafrostModel::dsaturations_dpc_liq(double s_i, double pc_liq,
        double pc_ice, double guess, double (&dsats)[3]) {

  if (!dsaturations_dpc_liq_if_above_freezing_(pc_liq, pc_ice, dsats)) {
    // temperature below 0
    DSatIce_DPClg_Functor_ func(s_i, pc_liq, pc_ice, wrm_);
    Tol_ tol(eps_);
    uintmax_t max_it(max_it_);
    double factor = 2.0;

    if (guess <= 0.) guess = 1.e-5;
    std::pair<double,double> result =
        boost::math::tools::bracket_and_solve_root(func, guess, factor, true, tol, max_it);

    ASSERT(max_it < max_it_);
    dsats[2] = result.first;

    // now liq and gas from ice
    dsats[1] = (1.0 - s_i)*wrm_->d_saturation(pc_liq) - dsats[2]*wrm_->saturation(pc_liq);
    dsats[0] = - dsats[1] - dsats[2];
  }
};

void WRMImplicitPermafrostModel::dsaturations_dpc_ice(double s_i, double pc_liq,
        double pc_ice, double guess, double (&dsats)[3]) {
  if (!dsaturations_dpc_ice_if_above_freezing_(pc_liq, pc_ice, dsats)) {
    // temperature below 0
    DSatIce_DPCil_Functor_ func(s_i, pc_liq, pc_ice, wrm_);
    Tol_ tol(eps_);
    uintmax_t max_it(max_it_);
    double factor = 2.0;

    if (guess <= 0.) guess = 1.e-5;
    std::pair<double,double> result =
        boost::math::tools::bracket_and_solve_root(func, guess, factor, true, tol, max_it);

    ASSERT(max_it < max_it_);
    dsats[2] = result.first;

    // now liq and gas from ice
    dsats[1] = - dsats[2]*wrm_->saturation(pc_liq);
    dsats[0] = - dsats[1] - dsats[2];
  }
};




} // namespace
} // namespace
} // namespace
