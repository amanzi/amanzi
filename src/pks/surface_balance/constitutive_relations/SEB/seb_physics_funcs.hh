/*
  Functions for calculating the snow / surface energy balance.
*/

#ifndef SURFACEBALANCE_SEB_PHYSICS_FUNCS_HH_
#define SURFACEBALANCE_SEB_PHYSICS_FUNCS_HH_

#include <cmath>
#include <string>
#include "seb_physics_defs.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace SEBPhysics {

// Main SEB functions
void UpdateIncomingRadiation(const SEB& seb, EnergyBalance& eb, bool debug=false);
void UpdateEvapResistance(const SEB& seb, EnergyBalance& eb, bool debug=false);
void UpdateEnergyBalance(const SEB& seb, const ThermoProperties& vp_surf,
                         EnergyBalance& eb, bool debug=false);

double DetermineSnowTemperature(const SEB& seb, ThermoProperties& vp_snow,
        EnergyBalance& eb, std::string method="toms");
void UpdateMassBalance(const SEB& seb, MassBalance& mb, EnergyBalance& eb,
                       SnowProperties& snow_new, bool debug=false);
void CalculateSurfaceBalance(SEB& seb, bool debug=false);

// Random helper functions
double CalcAlbedoSnow(double density_snow);
double CalcRoughnessFactor(double air_temp);

// Calculation of a snow temperature requires a root-finding operation, for
// which we use a functor.
class SnowTemperatureFunctor_ {
 public:

  explicit SnowTemperatureFunctor_(const SEB* seb,
          ThermoProperties* vp_snow,
          EnergyBalance* eb) :
      seb_(seb),
      vp_snow_(vp_snow),
      eb_(eb) {}

  double operator()(double temp) {
    vp_snow_->temp = temp;
    vp_snow_->UpdateVaporPressure();
    UpdateEnergyBalance(*seb_, *vp_snow_, *eb_);
    eb_->BalanceViaMelt();
    return eb_->fQm;
  }

 private:
  const SEB* const seb_;
  ThermoProperties* const vp_snow_;
  EnergyBalance* const eb_;
};


// Convergence criteria for root-finding
struct Tol_ {
  Tol_(double eps) : eps_(eps) {}
  bool operator()(const double& a, const double& b) const {
    return std::abs(a - b) <= eps_;
  }
  double eps_;
};

} // namespace
} // namespace
} // namespace

#endif
