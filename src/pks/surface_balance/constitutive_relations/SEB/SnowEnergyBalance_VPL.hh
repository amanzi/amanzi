/*
  Functions for calculating the snow-surface energy balance.

  Incoming Longwave radation is cacualted in this version, but if data is
  available we could incorporate it with the available met data.

  Atmospheric pressure is often used in snow models, If data is available we
  could incorperate it but for now Pa is held constant at 100 Pa.

  *** Equation for saturated vapor pressure over water is taken from Bolton,
      1980 'Monthly Weather Review'

  *** Equation for saturated vaport pressure over snow is taken from Buck,
      1996 'Buck Research Manual'

  *** See: http://cires.colorado.edu/~voemel/vp.html

*/

#ifndef SNOW_ENERGY_BALANCE_VPL_HH_
#define SNOW_ENERGY_BALANCE_VPL_HH_

namespace SurfaceEnergyBalance_VPL {

struct VaporPressure {
  double temp;
  double relative_humidity;
  double saturated_vaporpressure;
  double actual_vaporpressure;
  double dewpoint_temp;
};

struct EnergyBalance {
  double fQswIn;

  double Ps;                    // precip snow
  double Pr;                    // precip rain

  double temp_ground;           // ground temperature
  double temp_air;              // air temperature
  double temp_snow;
  double Us;                    // wind speed

  double dt;

  double water_depth;
  double water_fraction;
  double ht_snow;
  double density_snow;
  double age_snow;

  double air_vaporpressure;
  double snow_vaporpressure;
  double dewpoint_temp;
  double albedo_value;

  double stephB;
  double Apa;
  double SEs;
  double SEtun;
  double Dhe;
  double gZr;
  double rowaCp;
  double rowaLs;
  double rowaLe;
  double porrowaLe;
  double density_w;
  double density_freshsnow;
  double density_frost;
  //    double density_air;
  double Hf;
  double Ls;
  double Le;
  double VKc;
  //    double Cp;
  double Zr;
  double Zo;

  double fQlwIn;
  double QswIn;
  double fQlwOut;
  double fQh;
  double fQe;
  double fQc;
  double Qm;
  double Trw;
  double surface_pressure;
  double saturation_liquid;  
  double stored_surface_pressure;
  //double stored_fQe;
  double surface_porosity;
  double SurfaceVaporFlux;  // Second mass flux for cell center  
  double SWE;

  double Mr;
  double MIr;

  double AlbedoTrans;

};

struct LocalData {
  LocalData() {
    st_energy.stephB = 0.00000005670373;// Stephan-boltzmann Constant ------- [W/m^2 K^4]
    st_energy.Hf = 333500.0;            // Heat of fusion for melting snow -- [J/kg]
    st_energy.Ls = 2834000.0;           // Latent heat of sublimation ------- [J/kg]
    st_energy.Le = 2497848.;            // Latent heat of vaporization ------ [J/kg]
    st_energy.SEs = 0.98;               // Surface Emissivity for snow  ----- [-] ** From P. ReVelle (Thesis)
    st_energy.SEtun = 0.92;             // Surface Emissivity for tundra --- [-] ** From P. ReVelle (Thesis); Ling & Zhang, 2004
    st_energy.Zr = 10.0;                 // Referance ht of wind speed ------- [m]
 //   st_energy.Zo = 0.005;               // Roughness length  ---------------- [m] Mud flats, snow; no vegetation, no obstacles 
  //*Note on Roughness lenght* Should add Change from Snow 0.005 to bare ground 0.03 --> Open flat terrain; grass, few isolated obstacles.   
    st_energy.VKc = 0.41;               // Von Karman Constant -------------- [-]
    double Cp = 1004.0;                 // Specific heat of air ------------- [J/K kg]
    st_energy.Apa = 101.325;            // Atmospheric Pressure ------------- [KPa]

    st_energy.density_w = 1000;         // Density of Water ----------------- [kg/m^3]
    double density_air = 1.275;       // Density of Air ------------------- [kg/m^3]
    st_energy.density_frost = 200;      // Density of Frost (condensation) -- [kg/m^3]
    st_energy.density_freshsnow = 100;  // Density of Freshly fallebn snow -- [kg/m^3]

//    st_energy.density_frost = 250;
//    st_energy.density_freshsnow = 250;

    st_energy.gZr = 9.807*st_energy.Zr;
    st_energy.rowaCp = density_air*Cp;
    st_energy.rowaLs = density_air*st_energy.Ls;
    st_energy.rowaLe = density_air*st_energy.Le;

    vp_snow.relative_humidity = 1.;
    vp_ground.relative_humidity = 1.;
  }


  VaporPressure vp_air;
  VaporPressure vp_ground;
  VaporPressure vp_snow;
  EnergyBalance st_energy;

};


void UpdateIncomingRadiation(LocalData& seb);
void UpdateIncomingRadiationDerivatives(LocalData& seb);
void UpdateEFluxesSnow(LocalData& seb, double T);
double CalcMeltEnergy(LocalData& seb);
void UpdateGroundEnergy(LocalData& seb);
void UpdateGroundEnergyDerivatives(LocalData& seb);

void UpdateVaporPressure(VaporPressure& vp);
double CalcAlbedo(EnergyBalance& eb);
double EnergyBalanceResidual(LocalData& seb, double Xx);
double CalcSnowTemperature(LocalData& seb);

void UpdateMassMelt(EnergyBalance& eb);
void UpdateMassSublCond(EnergyBalance& eb);
void UpdateMassEvap(EnergyBalance& eb);
void WaterMassCorrection(EnergyBalance& eb);
void UpdateSnow(EnergyBalance& eb);

// Main "public" methods.
void SnowEnergyBalance(LocalData& seb);
void UpdateEnergyBalance(LocalData& seb);
void UpdateEnergyBalanceDerivative(LocalData& seb);

}// Namespace

#endif
