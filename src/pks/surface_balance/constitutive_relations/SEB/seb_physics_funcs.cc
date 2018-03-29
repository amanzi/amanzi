/*
  Functions for calculating the snow / surface energy balance.
*/

#include <iostream>
#include <cmath>
#include <algorithm>
#include "boost/math/tools/roots.hpp"

#include "dbc.hh"

#include "seb_physics_funcs.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace SEBPhysics {

#define SWE_EPS 1.e-12
#define ENERGY_BALANCE_TOL 1.e-8


double CalcAlbedoSnow(double density_snow) {
  double AlSnow;
  if (density_snow <= 432.23309912785146) {
    AlSnow = 1.0 - 0.247 * std::pow(0.16 + 110*std::pow(density_snow/1000, 4), 0.5);
  } else {
    AlSnow = 0.6 - density_snow / 4600;
  }
  return AlSnow;
}

double CalcRoughnessFactor(double snow_height, double Z_rough_bare, double Z_rough_snow)
{
  double Zfraction = snow_height <= 0. ? 0. : snow_height >= Z_rough_bare ? 1. :
                     1 - (Z_rough_bare-snow_height)/Z_rough_bare;
  return Z_rough_snow * Zfraction + Z_rough_bare * (1-Zfraction);
}


std::pair<double,double> IncomingRadiation(const MetData& met, double albedo)
{
  // Calculate incoming short-wave radiation
  double fQswIn = (1 - albedo) * met.QswIn;

  // Calculate incoming long-wave radiation
  double fQlwIn = met.QlwIn;

  return std::make_pair(fQswIn, fQlwIn);
}


//
// Calculate longwave from air temp and relative humidity
// ------------------------------------------------------------------------------------------
double CalcIncomingLongwave(double air_temp, double relative_humidity, double c_stephan_boltzmann) {
  double e_air = std::pow(10 * VaporPressureAir(air_temp, relative_humidity), air_temp / 2016.);
  e_air = 1.08 * (1 - std::exp(-e_air));
  return e_air * c_stephan_boltzmann * std::pow(air_temp,4);
}


double OutgoingRadiation(double temp, double emissivity, double c_stephan_boltzmann)
{
  // Calculate outgoing long-wave radiation
  return emissivity * c_stephan_boltzmann * std::pow(temp,4);
}

double WindFactor(double Us, double Z_Us, double Z_rough, double c_von_Karman)
{
  // Calculate D_h, D_e, 
  return std::pow(c_von_Karman,2) * Us / std::pow(std::log(Z_Us / Z_rough), 2);
}


double StabilityFunction(double air_temp, double skin_temp, double Us,
                         double Z_Us, double c_gravity)
{
  double Ri  = c_gravity * Z_Us * (air_temp - skin_temp) / (air_temp * std::pow(Us,2));
  if (Ri >= 0.) {
    // stable condition
    return 1. / (1 + 10*Ri);
  } else {
    // Unstable condition
    return (1 - 10*Ri);
  }
}


double SaturatedVaporPressure(double temp)
{
  // Sat vap. press o/water Dingman D-7 (Bolton, 1980)
  // *** (Bolton, 1980) Calculates vapor pressure in [kPa]  ****
  double tempC = temp - 273.15;
  return 0.6112 * std::exp(17.67 * tempC / (tempC + 243.5));
}

double VaporPressureAir(double air_temp, double relative_humidity)
{
  return SaturatedVaporPressure(air_temp) * relative_humidity;
}

double VaporPressureGround(const GroundProperties& surf, const ModelParams& params)
{
  // Ho & Webb 2006
  double relative_humidity = -1;
  if (surf.pressure < params.Apa * 1000.) {
    // vapor pressure lowering
    double pc = 1000.*params.Apa - surf.pressure;
    relative_humidity = std::exp(-pc / (surf.density_w * params.R_ideal_gas * surf.temp));
  } else {
    relative_humidity = 1.;
  }
  return relative_humidity * SaturatedVaporPressure(surf.temp);
}


double EvaporativeResistanceGround(const GroundProperties& surf,
        const MetData& met,
        const ModelParams& params, 
        double vapor_pressure_air, double vapor_pressure_ground)
{
  // calculate evaporation prefactors
  if (vapor_pressure_air > vapor_pressure_ground) { // condensation
    return 0.;
  } else if (surf.pressure > 1000.*params.Apa) { // ponded water present
    return 0.;
  } else { // evaporating from soil
    // Equation for reduced vapor diffusivity
    // See Sakagucki and Zeng 2009 eqaution (9) and Moldrup et al., 2004. 
    double vp_diffusion = 0.000022 * (std::pow(surf.porosity,2))
                          * std::pow((1-(0.0556/surf.porosity)),(2+3*params.Clapp_Horn_b));

    // Sakagucki and Zeng 2009 eqaution (10)
    double L_Rsoil = std::exp(std::pow(surf.saturation_gas, 5));
    L_Rsoil = surf.dz * (L_Rsoil -1) * (1/(std::exp(1.)-1));
    double Rsoil = L_Rsoil/vp_diffusion;
    return Rsoil;
  }
}


double SensibleHeat(double resistance_coef,
                    double density_air,
                    double Cp_air,
                    double air_temp,
                    double skin_temp) {
  return resistance_coef * density_air * Cp_air * (air_temp - skin_temp);
}


double LatentHeat(double resistance_coef,
                  double density_air,  /// this should be w?
                  double latent_heat_fusion,
                  double vapor_pressure_air,
                  double vapor_pressure_skin,
                  double Apa) {
  return resistance_coef * density_air * latent_heat_fusion * 0.622
      * (vapor_pressure_air - vapor_pressure_skin) / Apa;
}

double ConductedHeatIfSnow(double ground_temp,
                           const SnowProperties& snow)
{
  // Calculate heat conducted to ground, if snow
  double Ks = -1;
  if (snow.density > 150) { // frost hoar
    double snow_hoar_density = 1. / ((0.90/snow.density) + (0.10/150));
    Ks = 2.9e-6 * std::pow(snow_hoar_density,2);
  } else {
    Ks = 2.9e-6 * std::pow(snow.density,2);
  }
  return Ks * (snow.temp - ground_temp) / snow.height;
}


void UpdateEnergyBalanceWithSnow(const GroundProperties& surf,
        const SnowProperties& snow,
        const MetData& met,
        const ModelParams& params,
        EnergyBalance& eb)
{
  // incoming radiation -- DONE IN MAIN

  // outgoing radiation
  eb.fQlwOut = OutgoingRadiation(snow.temp, snow.emissivity, params.stephB);

  // sensible heat
  double Dhe = WindFactor(met.Us, met.Z_Us, CalcRoughnessFactor(snow.height, surf.roughness, snow.roughness), params.VKc);
  double Sqig = StabilityFunction(met.air_temp, snow.temp, met.Us, met.Z_Us, params.gravity);
  eb.fQh = SensibleHeat(Dhe * Sqig, params.density_air, params.Cp_air, met.air_temp, snow.temp);

  // latent heat
  double vapor_pressure_air = VaporPressureAir(met.air_temp, met.relative_humidity);
  double vapor_pressure_skin = SaturatedVaporPressure(snow.temp);
  eb.fQe = LatentHeat(Dhe * Sqig, params.density_air, params.Ls, vapor_pressure_air, vapor_pressure_skin,
                      params.Apa);

  // conducted heat
  eb.fQc = ConductedHeatIfSnow(surf.temp, snow);

  // balance of energy goes into melting
  eb.fQm = eb.fQswIn + eb.fQlwIn - eb.fQlwOut + eb.fQh - eb.fQc;
  // only add latent heat if evaporating (condensation's energy gets lost to atmosphere, not surface)
  if (eb.fQe < 0) eb.fQm += eb.fQe;
}


void UpdateEnergyBalanceWithoutSnow(const GroundProperties& surf,
        const MetData& met,
        const ModelParams& params,
        EnergyBalance& eb)
{
  // incoming radiation
  std::tie(eb.fQswIn, eb.fQlwIn) = IncomingRadiation(met, surf.albedo);

  // outgoing radiation
  eb.fQlwOut = OutgoingRadiation(surf.temp, surf.emissivity, params.stephB);

  // sensible heat
  double Dhe = WindFactor(met.Us, met.Z_Us, surf.roughness, params.VKc);
  double Sqig = StabilityFunction(met.air_temp, surf.temp, met.Us, met.Z_Us, params.gravity);
  eb.fQh = SensibleHeat(Dhe*Sqig, params.density_air, params.Cp_air, met.air_temp, surf.temp);

  // latent heat
  double vapor_pressure_air = VaporPressureAir(met.air_temp, met.relative_humidity);
  double vapor_pressure_skin = VaporPressureGround(surf, params);

  double Rsoil = EvaporativeResistanceGround(surf, met, params, vapor_pressure_air, vapor_pressure_skin);
  double coef = 1.0 / (Rsoil + 1.0/(Dhe*Sqig));
  if (surf.temp < 273.15) {
    eb.fQe = LatentHeat(coef, params.density_air, params.Ls,
                        vapor_pressure_air, vapor_pressure_skin, params.Apa);
  } else {
    eb.fQe = LatentHeat(coef, params.density_air, params.Le,
                        vapor_pressure_air, vapor_pressure_skin, params.Apa);
  }
  
  // melt, conducted is 0
  eb.fQm = 0.;
  eb.fQc = eb.fQswIn + eb.fQlwIn - eb.fQlwOut + eb.fQh;
}

// Snow temperature calculation.
double DetermineSnowTemperature(const GroundProperties& surf,
        const MetData& met,
        const ModelParams& params, 
        SnowProperties& snow,
        EnergyBalance& eb,
        std::string method)
{
  SnowTemperatureFunctor_ func(&surf, &snow, &met, &params, &eb);
  Tol_ tol(ENERGY_BALANCE_TOL);
  boost::uintmax_t max_it(100);
  double left, right;
  double res_left, res_right;

  double res_init = func(surf.temp);
  if (res_init < 0.) {
    right = surf.temp;
    res_right = res_init;

    left = surf.temp - 1.;
    res_left = func(left);
    while (res_left < 0.) {
      right = left;
      res_right = res_left;
      left = left - 1.;
      res_left = func(left);
    }
  } else {
    left = surf.temp;
    res_left = res_init;

    right = surf.temp + 1.;
    res_right = func(right);
    while (res_right > 0.) {
      left = right;
      res_left = res_right;
      right = right + 1.;
      res_right = func(right);
    }
  }

  //  std::cout << "Determining snow temp in interval (" << left << "," << right << ")" << std::endl;
  
  std::pair<double,double> result;
  auto my_max_it = max_it;
  if (method == "bisection") {
    result = boost::math::tools::bisect(func, left, right, tol, max_it);
  } else if (method == "toms") {
    result = boost::math::tools::toms748_solve(func, left, right, res_left, res_right, tol, max_it);
  }

  if (max_it >= my_max_it) throw("Nonconverged Surface Energy Balance");
  // boost algorithms calculate the root such that the root is contained within the interval
  // [result.first, result.second], and that width of that interval is less than TOL.
  // We choose to set the solution as the center of that interval.
  // Call the function again to set the fluxes.
  double solution = (result.first + result.second)/2.;
  //  std::cout << "  Got T=" << solution << " with SW = " << eb.fQswIn << ", LW = " << eb.fQlwIn << ", LWout = " << eb.fQlwOut << ", Sens = " << eb.fQh << ", Lat = " << eb.fQe << ", Cond = " << eb.fQc << ", melt = " << eb.fQm << std::endl;

  return solution;
}



MassBalance UpdateMassBalance(const GroundProperties& surf, const SnowProperties& snow_old,
        const MetData& met, const ModelParams& params, const EnergyBalance eb,
        SnowProperties& snow_new, double dt)
{
  MassBalance mb;
  mb.dt = dt;

  // Melt rate given by available energy rate divided by heat of fusion.
  mb.Mm = eb.fQm / (surf.density_w * params.Hf);

  // Snow balance
  if (snow_old.height > 0.) {
    mb.Me = eb.fQe / (surf.density_w * params.Ls);

    double swe_old = snow_old.height * snow_old.density / surf.density_w;
    double swe_new = swe_old + (met.Ps - mb.Mm + mb.Me)*mb.dt;
    // First do a pass to ensure we are not melting or sublimating ALL
    // of the available snow.  If so, adjust dt
    if (swe_new < 0.) {
      // Must adjust... we are melting or sublimating all of the available snow.
      ASSERT(met.Ps - mb.Mm + mb.Me < 0.);
      mb.dt = -swe_old / (met.Ps - mb.Mm + mb.Me);
      swe_new = 0.;
    }

    // All future calculations work with the new dt, which will
    // exactly result in 0 snow height (if it would have been negative).

    // age the old snow
    double age_settled = snow_old.age + mb.dt / 86400.;
    double dens_settled = params.density_freshsnow
        * std::max(std::pow(age_settled, 0.3), 1.);

    // Match frost age with assigned density -- Calculate which day frost
    // density matched snow defermation function from (Martinec, 1977)
    double age_frost = std::pow((params.density_frost / params.density_freshsnow),
				(1/0.3)) - 1 + mb.dt / 86400.;

    // precip
    double age_precip = mb.dt / 86400.;

    // determine the new height
    // -- sources
    double swe_settled = swe_old;
    double swe_frost = mb.Me > 0. ? mb.Me*mb.dt : 0.;
    double swe_precip = met.Ps*mb.dt;

    // -- sinks
    double swe_subl =  mb.Me < 0. ? -mb.Me*mb.dt : 0.;
    double swe_melt = mb.Mm * mb.dt;

    // -- sublimate precip first
    ASSERT(swe_subl >= 0.);
    if (swe_subl > 0.) {
      if (swe_subl > swe_precip) {
	swe_subl -= swe_precip;
	swe_precip = 0.;
      } else {
	swe_precip -= swe_subl;
	swe_subl = 0.;
      }
    }

    // -- next sublimate settled snow
    if (swe_subl > 0.) {
      ASSERT(swe_subl <= swe_settled + SWE_EPS);
      swe_settled -= swe_subl;
      swe_subl = 0.;
    }

    // -- melt settled snow first
    ASSERT(swe_melt >= -SWE_EPS);
    if (swe_melt > 0.) {
      if (swe_melt > swe_settled) {
        swe_melt -= swe_settled;
        swe_settled = 0.;
      } else {
        swe_settled -= swe_melt;
        swe_melt = 0.;
      }
    }

    // -- now melt frost, precip by even amounts
    if (swe_melt > SWE_EPS) {
      ASSERT(swe_frost + swe_precip > 0.);
      double swe_melt_from_frost = swe_melt * (swe_frost / (swe_frost + swe_precip));
      double swe_melt_from_precip = swe_melt - swe_melt_from_frost;

      swe_frost -= swe_melt_from_frost;
      swe_precip -= swe_melt_from_precip;
    }

    // -- check we didn't screw up
    ASSERT(swe_settled >= -SWE_EPS);
    ASSERT(swe_frost >= -SWE_EPS);
    ASSERT(swe_precip >= -SWE_EPS);

    // -- convert these to heights
    double ht_settled = swe_settled * surf.density_w / dens_settled;
    double ht_frost = swe_frost * surf.density_w / params.density_frost;
    double ht_precip = swe_precip * surf.density_w / params.density_freshsnow;

    // set the snow properties
    double swe_total = std::max(swe_settled + swe_frost + swe_precip, 0.);
    ASSERT(std::abs(swe_total - swe_new) < SWE_EPS);
    snow_new.height = std::max(ht_settled + ht_frost + ht_precip, 0.);
    snow_new.age = swe_new > 0. ? (swe_settled*age_settled + swe_frost*age_frost + swe_precip*age_precip) / swe_new : 0.;
    snow_new.density = snow_new.height > 0. ? swe_new * surf.density_w / snow_new.height : params.density_freshsnow;
    snow_new.SWE = swe_new;

    // set the water properties
    // -- water source to ground is (corrected) melt and rainfall
    // NOTE: these rates can only be correct if over mb.dt
    

  } else {
    // set the evaporative flux of mass
    if (surf.temp < 273.15) {
      mb.Me = eb.fQe / (surf.density_w * params.Ls);
    } else {
      mb.Me = eb.fQe / (surf.density_w * params.Le);
    }      

    // set the snow properties
    snow_new.height = met.Ps * mb.dt
        * surf.density_w / params.density_freshsnow;
    snow_new.age = mb.dt / 86400.;
    snow_new.density = params.density_freshsnow;
    snow_new.SWE = snow_new.height * snow_new.density / surf.density_w;
  }

  // if (debug) {
  //   std::cout << "Mass Balance:\n"
  //             << "  Mm   = " << mb.Mm << std::endl
  //             << "  Me   = " << mb.Me << std::endl
  //             << "  Ps   = " << met.Ps << std::endl
  //             << "  Pr   = " << met.Pr << std::endl
  //             << "  Snow Melt:\n"
  //             << "    old ht   = " << snow_old.height << std::endl
  //             << "    new ht   = " << snow_new.height << std::endl
  //             << "    new age  = " << snow_new.age << std::endl
  //             << "    new dens = " << snow_new.density << std::endl
  //             << "    SWE      = " << snow_new.SWE << std::endl
  //             << "  Water Balance:\n"
  //             << "    surf src = " << mb.MWg << std::endl
  //             << "    sub src  = " << mb.MWg_subsurf << std::endl;
  // }
  return mb;
}



// master driver
std::tuple<SnowProperties, EnergyBalance, MassBalance, FluxBalance>
CalculateSurfaceBalance(double dt,
                        const GroundProperties& surf,
                        const SnowProperties& snow_old,
                        const MetData& met,
                        const ModelParams& params,
                        bool debug, const Teuchos::RCP<VerboseObject>& vo)
{
  EnergyBalance eb;
  SnowProperties snow_new(snow_old);

  if (snow_new.height > 0.) {
    // snow on the ground, solve for snow temperature
    snow_new.albedo = CalcAlbedoSnow(snow_new.density);
    std::tie(eb.fQswIn, eb.fQlwIn) = IncomingRadiation(met, snow_new.albedo);
    snow_new.temp = DetermineSnowTemperature(surf, met, params, snow_new, eb);

    if (snow_new.temp > 273.15) {
      // limit snow temp to 0, then melt with the remaining energy
      snow_new.temp = 273.15;
      UpdateEnergyBalanceWithSnow(surf, snow_new, met, params, eb);
      eb.error = 0.;
      
    } else {
      // snow not melting
      UpdateEnergyBalanceWithSnow(surf, snow_new, met, params, eb);
      eb.error = eb.fQm;
      eb.fQm = 0.;
    }

  } else {
    // no snow on the ground, balance given by conduction
    UpdateEnergyBalanceWithoutSnow(surf, met, params, eb);
  }

  // Mass balance
  MassBalance mb = UpdateMassBalance(surf, snow_old, met, params, eb, snow_new, dt);

  // update fluxes delivered to ground, except evaporation
  FluxBalance flux;
  // mass to surface is melt and precip, with melt normalized to dt
  flux.M_surf = met.Pr + mb.Mm * mb.dt / dt;

  // energy to surface is conducted energy, enthalpy of rain and melt
  flux.E_surf = eb.fQc + surf.density_w * met.Pr * (met.air_temp-273.15) * params.Cv_water;
  //           + surf.density_w * met.Mm * 0. * params.Cv_water; // enthalpy of water at 0c is 0

  // mass, energy to subsurface is 0
  flux.M_sub = 0.;
  flux.E_sub = 0.;

  // now add in evaporation  
  if (snow_new.height == 0.) {
    // check to see if evaporation takes all the water
    double evap_to_surface = 1.;
    if (mb.Me > 0 && surf.ponded_depth < params.evap_transition_width) {
      evap_to_surface = surf.ponded_depth / params.evap_transition_width;
    }
    flux.M_surf += evap_to_surface * mb.Me;
    flux.M_sub += (1-evap_to_surface) * mb.Me;

    if (eb.fQe > 0.) {
      // if condensation, just add energy due to introduction of new water.
      // Latent heat released goes into the atmosphere.
      flux.E_surf += evap_to_surface * mb.Me * surf.density_w * (met.air_temp-273.15) * params.Cv_water;
      flux.E_sub += (1-evap_to_surface) * mb.Me * surf.density_w * met.air_temp * params.Cv_water;
    } else {
      // if evaporation, take away the energy of that water, but also cool due
      // to latent heat removal.  Note this is equivalent to simply taking
      // away the energy of the VAPOR removed.
      flux.E_surf += evap_to_surface * (mb.Me * surf.density_w * (surf.temp-273.15) * params.Cv_water + eb.fQe);
      flux.E_sub += (1-evap_to_surface) * (mb.Me * surf.density_w * surf.temp * params.Cv_water + eb.fQe);
    }
    
  }
  
  if (vo.get()) {
    *vo->os() << "---------------------------------------------------------" << std::endl
              << "Surface Energy Balance:" << std::endl
              << "  Incoming Radiation Energy Terms:" << std::endl
              << "    windspeed, Z: " << met.Us << "  " << met.Z_Us << std::endl
              << "    fQswIn   = " << eb.fQswIn << std::endl
              << "    fQlwIn   = " << eb.fQlwIn << std::endl
              << "  Evap/Cond Terms:" << std::endl
              << "    air temp, skin temp: " << met.air_temp << "  " << (snow_old.height > 0 ? snow_new.temp : surf.temp) << std::endl
              << "    skin saturation_gas: " << surf.saturation_gas << std::endl
              << "  Energy Balance Terms (ht_snow = " << snow_old.height << "):" << std::endl
              << "    SnowSurfaceTemp  = " << snow_new.temp << std::endl
              << "    GroundSurfaceTemp  = " << surf.temp << std::endl
              << "    fQlwOut  = " << eb.fQlwOut << std::endl
              << "    fQh      = " << eb.fQh << std::endl
              << "    fQe      = " << eb.fQe << std::endl
              << "    fQc      = " << eb.fQc << std::endl
              << "    fQm      = " << eb.fQm << std::endl
              << "    eb error = " << eb.error << std::endl
              << "  Mass Balance:\n"
              << "    Mm   = " << mb.Mm << std::endl
              << "    Me   = " << mb.Me << std::endl
              << "    Ps   = " << met.Ps << std::endl
              << "    Pr   = " << met.Pr << std::endl
              << "    Snow Melt:\n"
              << "      old ht   = " << snow_old.height << std::endl
              << "      new ht   = " << snow_new.height << std::endl
              << "      old age  = " << snow_old.age << std::endl
	      << "      new age  = " << snow_new.age << std::endl
	      << "      old dens = " << snow_old.density << std::endl
              << "      new dens = " << snow_new.density << std::endl
              << "      SWE      = " << snow_new.SWE << std::endl
              << "  Fluxes to other layeres:\n"
              << "      surf mass src = " << flux.M_surf << std::endl
              << "      sub mass src = " << flux.M_sub << std::endl
              << "      surf energy src = " << flux.E_surf << std::endl
              << "      sub energysrc = " << flux.E_sub << std::endl;
  }
  return std::make_tuple(snow_new, eb, mb, flux);
}


} // namespace
} // namespace
} // namespace
