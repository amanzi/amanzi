/*
  Data structures and parameters for calculating the snow / surface energy balance.

  Nomenclature: A few names are used which imply various surface processes.
    -- snow, water, ice, and tundra are the four possible surface materials
    -- pond indicates ponded water, either ice or water
    -- ground indicates whatever is below any snow -- ice or water or tundra
*/


#ifndef SURFACEBALANCE_SEB_PHYSICS_DEFS_HH_
#define SURFACEBALANCE_SEB_PHYSICS_DEFS_HH_

#include <limits>

namespace Amanzi {
namespace SurfaceBalance {
namespace SEBPhysics {


// Struct of thermodynamic properties
struct ThermoProperties {
  double temp;                          // temperature [K]
  double pressure;                      // [Pa]
  double porosity;                      // [-]
  double density_w;                     // density [kg/m^3]
  double relative_humidity;             // [-]
  double saturated_vaporpressure;       // [kPa]
  double actual_vaporpressure;          // [kPa]
  double dewpoint_temp;                 // [K]

  ThermoProperties() :
      temp(std::numeric_limits<double>::signaling_NaN()),
      pressure(101325.),
      porosity(1.),
      density_w(1000.),
      relative_humidity(std::numeric_limits<double>::signaling_NaN()),
      saturated_vaporpressure(std::numeric_limits<double>::signaling_NaN()),
      actual_vaporpressure(std::numeric_limits<double>::signaling_NaN()),
      dewpoint_temp(std::numeric_limits<double>::signaling_NaN()) {}

  void UpdateVaporPressure();
};


// Struct of snow state
struct SnowProperties {
  double ht;                    // snow depth [m] (NOT SWE!)
  double density;               // snow density [ kg / m^3 ]
  double age;                   // snow age [days]

  SnowProperties() :
      ht(std::numeric_limits<double>::signaling_NaN()),
      density(std::numeric_limits<double>::signaling_NaN()),
      age(std::numeric_limits<double>::signaling_NaN()) {}
};


// Struct of surface properties for whatever is exposed to air.
struct SurfaceProperties {
  double albedo;
  double Zo;
  double emissivity;

  SurfaceProperties() :
      albedo(std::numeric_limits<double>::signaling_NaN()),
      Zo(std::numeric_limits<double>::signaling_NaN()),
      emissivity(std::numeric_limits<double>::signaling_NaN()) {}
};


// struct of input MetData.
struct MetData {
  double Us;                    // wind speed, [m/s]
  double QswIn;                 // incoming short-wave radiation, [J/ ??]
  double Ps;                    // precip snow, [m (SWE)/s]
  double Pr;                    // precip rain, [m/s]
  ThermoProperties vp_air;

  MetData() :
      Us(std::numeric_limits<double>::signaling_NaN()),
      QswIn(std::numeric_limits<double>::signaling_NaN()),
      Ps(std::numeric_limits<double>::signaling_NaN()),
      Pr(std::numeric_limits<double>::signaling_NaN()),
      vp_air() {}
};


// Catch-all of leftover parameters.
struct ModelParams {
  // densities of water in various forms
  double density_air;
  double density_freshsnow, density_frost;

  // thermodynamic constants of water
  double Hf, Ls, Le, Cp;

  // constants for energy equations
  double VKc;
  double Zr;
  double stephB;

  // other constants
  double Apa;
  double evap_transition_width;
  double gravity;

  ModelParams() :
      density_air(1.275),       // [kg/m^3]
      density_freshsnow(100.),  // [kg/m^3]
      density_frost(200.),      // [kg/m^3]
      Hf(333500.0),             // Heat of fusion for melting snow -- [J/kg]
      Ls(2834000.0),            // Latent heat of sublimation ------- [J/kg]
      Le(2497848.),             // Latent heat of vaporization ------ [J/kg]
      Cp(1004.0),               // Specific heat of air ------------- [J/K kg]
      VKc(0.41),                // Von Karman Constant -------------- [-]
      Zr(2.0),                  // Reference height of wind speed --- [m]
      stephB(0.00000005670373), // Stephan-Boltzmann constant ------- [W/m^2 K^4]
      Apa(101.325),             // atmospheric pressure ------------- [kPa]
      evap_transition_width(100.), // transition on evaporation from surface to evaporation from subsurface [m]
      gravity(9.807) {}         // gravity [kg m / s^2]
};


// Global ModelInput struct -- all input to SEB
struct ModelInput {
  double dt;                    // time step size [s]

  ThermoProperties vp_ground;      // vapor pressure properties of soil/ponded water
  ThermoProperties vp_snow;      // vapor pressure properties of snow
  SurfaceProperties surf;
  SnowProperties snow_old;
  MetData met;

  ModelInput() :
      dt(std::numeric_limits<double>::signaling_NaN()),
      vp_ground(),
      vp_snow(),
      surf(),
      snow_old(),
      met() {}
};


// Struct collecting energy balance terms.
struct EnergyBalance {  // all are [J/ (m^2 s)]
  double fQswIn;        // incoming short-wave radiation
  double fQlwIn;        // incoming long-wave radiation
  double fQlwOut;       // outgoing long-wave radiation
  double fQh;           // sensible heat
  double fQe;           // latent heat
  double fQc;           // heat conducted to ground surface
  double fQm;           // energy available for melting snow
  double Dhe;           // special constant for use in e and h, precalculated for efficiency

  EnergyBalance() :
      fQswIn(std::numeric_limits<double>::signaling_NaN()),
      fQlwIn(std::numeric_limits<double>::signaling_NaN()),
      fQlwOut(std::numeric_limits<double>::signaling_NaN()),
      fQh(std::numeric_limits<double>::signaling_NaN()),
      fQe(std::numeric_limits<double>::signaling_NaN()),
      fQc(std::numeric_limits<double>::signaling_NaN()),
      fQm(std::numeric_limits<double>::signaling_NaN()),
      Dhe(std::numeric_limits<double>::signaling_NaN()) {}

  void BalanceViaMelt();
  void BalanceViaConduction();
};


// Struct collecting mass balance terms.
struct MassBalance {    // all are in [m/s] of WATER, i.e. snow are in SWE
  double Me;    // condensation of water/frost (if positive),
                // sublimation/evaporation of snow/water (if negative)
  double Mm;    // melt rate (positive indicates increasing water, decreasing snow)
  double MWg;   // water source to ground
  double MWg_subsurf; // water source to subsurface cell
  double MWg_temp; // temperature of sources

  MassBalance() :
      Me(std::numeric_limits<double>::signaling_NaN()),
      Mm(std::numeric_limits<double>::signaling_NaN()),
      MWg(std::numeric_limits<double>::signaling_NaN()) {}
};


// Global struct of all model output.
struct ModelOutput {
  EnergyBalance eb;
  MassBalance mb;
  SnowProperties snow_new;
};


// Global struct of all data.
struct SEB {
  // model input -- these should never change
  ModelInput in;

  // model output -- required by ATS
  ModelOutput out;

  // parameters -- physical constants, these should never change
  ModelParams params;
};


// Used to calculate surface properties, prior to calling SEB.
struct SurfaceParams {
  double a_tundra, a_water, a_ice;      // albedos
  double e_snow, e_tundra, e_water, e_ice;     // emissivities
  double Zsmooth, Zrough;      // roughness coefs

  SurfaceParams() :
      a_tundra(0.15),           // [-] Grenfell and Perovich, (2004)
      a_water(0.141),           // [-] Cogley J.G. (1979)
      a_ice(0.44),              // [-] deteriorated ice, Grenfell and Perovich, (2004)
      e_snow(0.98),             // [-] emissivity for snow, From P. ReVelle (Thesis)
      e_tundra(0.92),           // [-] emissivity for tundra, From P. ReVelle
                                //         (Thesis), Ling & Zhang, 2004
      e_water(0.995),           // [-] emissivity of water >> tundra, this is from Wiki
      e_ice(0.98),              // [-] emissivity of ice >> tundra, this is from Wiki
      Zsmooth(0.005),           // [m]? roughness coef of smooth
      Zrough(0.03) {}           // [m]? roughness coef of rough
};


// Partitioning of weights for smoothing surface properties, output data structure.
struct Partition {
  double perSnow;
  double perWater;
  double perIce;
  double perTundra;

  double Interpolate(double snow, double water, double ice, double tundra) {
    return snow*perSnow + water*perWater + ice*perIce + tundra*perTundra;
  }
};


// Partitioning of weights for smoothing surface properties.
struct Partitioner {
  double water_pen, snow_pen;   // radiation penetration depths

  Partitioner() :
      water_pen(0.1),           // [m] QswIn penetration depth through water
      snow_pen(0.02) {}         // [m] QswIn penetration depth through snow
  Partitioner(double water_pen_, double snow_pen_) :
      water_pen(water_pen_),    // [m] QswIn penetration depth through water
      snow_pen(snow_pen_) {}    // [m] QswIn penetration depth through snow

  Partition CalcPartition(double ht_snow, double ht_pond,
                          double unfrozen_fraction) const;
};


} // namespace
} // namespace
} // namespace


#endif
