/*
  Data structures and parameters for calculating the snow / surface energy balance.
*/


#include <cmath>

#include "seb_physics_defs.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace SEBPhysics {

void ThermoProperties::UpdateVaporPressure() {  // Ho & Webb 2006
  double pressure_fudgefactor = 100;
  double R = 461.52; // Pa m^3 / kg K
  if (std::isnan(relative_humidity)) {
    if (pressure < 101325.) {
      // vapor pressure lowering
      double pc = 101325. - pressure;
      relative_humidity = std::exp(-pc*pressure_fudgefactor / (density_w*R*temp));
    } else {
      relative_humidity = 1.;
    }
  }

  double tempC = temp - 273.15;
  // Sat vap. press o/water Dingman D-7 (Bolton, 1980)
  // *** (Bolton, 1980) Calculates vapor pressure in [kPa]  ****
  saturated_vaporpressure = 0.6112*std::exp(17.67*tempC / (tempC+243.5));
  actual_vaporpressure = saturated_vaporpressure * relative_humidity;
}


Partition Partitioner::CalcPartition(double ht_snow,
        double ht_pond, double unfrozen_fraction) const {
  Partition part;

  if (ht_snow > snow_pen) {
    // incoming radiation does not penetrate the snow
    part.perSnow = 1.;
    part.perTundra = 0.;
    part.perWater = 0.;
    part.perIce = 0.;

  } else {
    part.perSnow = std::pow(ht_snow/snow_pen, 2);
    double remainder = 1. - part.perSnow;

    if (ht_pond > water_pen) {
      // incoming radiation does not penetrate snow and ponding
      part.perWater = unfrozen_fraction * remainder;
      part.perIce = remainder - part.perWater;
      part.perTundra = 0.;
    } else {
      double perPond = remainder * ht_pond / water_pen;
      part.perWater = unfrozen_fraction * perPond;
      part.perIce = (1. - unfrozen_fraction) * perPond;
      part.perTundra = remainder - perPond;
    }
  }
  return part;
}

void EnergyBalance::BalanceViaMelt() {
  // Melt energy is the balance
  fQm = fQswIn + fQlwIn + fQlwOut + fQh + fQe - fQc;
}

void EnergyBalance::BalanceViaConduction() {
  // Melt energy is the balance
  fQc = fQswIn + fQlwIn + fQlwOut + fQh + fQe;
}



} // namespace
} // namespace
} // namespace
