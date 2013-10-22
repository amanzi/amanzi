#ifndef SNOW_ENERGY_BALANCE_
#define SNOW_ENERGY_BALANCE_
#include <sstream>
#include "SnowEnergyBalance.hh"

/* THIS CODE WILL CALCULATE THE SNOWSURFACE ENERGY BALANCE THE SNOW SURFRACE TEMPERUATURE.

**** Incomming Longwave radation is cacualted in this version, but if data is avaialbe we could incoroerate it with the avialable met data****

+++Atmospheric pressure is often used in snow models, If data is available we could incorperate it but for now Pa is held constant at
100 pa++++

*** Equation for saturated vapor pressure over water is taken from Bolton, 1980 'Monthly Weather Review'
*** Equation for saturated vaport pressure over snow is taken from Buck, 1996 'Buck Research Manual'
See: http://cires.colorado.edu/~voemel/vp.html
***************************** */

namespace SurfaceEnergyBalance {

struct VaporPressure {
  double temp;
  double relative_humidity;
  double saturated_vaporpressure;
  double actual_vaporpressure;
  double dewpoint_temp;
};

struct EnergyBalance {
  double air_temp;
  double relative_rumidity;
  double Us;
  double fQswIn;
  double Ps;
  double Pr;

  double Tb;
  double Dt;
  double water_depth;
  double por;

  double ht_snow;
  double HoldZs;
  double density_snow;
  double HoldROWs;
  double nosnowdays;

  double air_vaporpressure;
  double snow_vaporpressure;
  double dewpoint_temp;
  double snow_conduct_value;
  double albedo_value;

  double stephB;
  double Apa;
  double SEs;
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
  double Ts;
  double Trw;

  double Ml;
  double Mr;
  double SublL;
  double SublR;
  double EvR;
  double EvL;
  double CiL;
  double TotwLoss;

  double varvar;
  std::string funcall;


};

//vapor_pressure * vp_pointer;
struct LocalData {
  VaporPressure vp_air;
  VaporPressure vp_ground;
  VaporPressure vp_snow;
  EnergyBalance st_energy;

};

void CalcEFluxTempIndependent (LocalData& seb);
void CalcEFluxTempDependent (LocalData& seb, double T);

// #### FUNCTIONS TO CALCULATE SATURATION VAPOR PRESSURE & ACTUALL VAPOR PRESSURE
void VaporCalc (VaporPressure& vp);


// #### FUNCTIONS TO CALCULATE ALBEDO
void AlbedoCalc (EnergyBalance& eb);


// ### FUNCTION TO CALCULATE THERMAL CONDUCTIVITY OF SNOW
void ThermalConductSnow (EnergyBalance& eb);


/*  FUNCTION TO SOLVE SNOW SURFACE ENERGY BALANCE FOR SURFACE TEMP (Ts)
    ###############################################################################################*
    Bisection Method
    (1-albedo)QswIn + Qlwin + QlwOut(Ts) + Qh(Ts) + Qc(Ts) +  Qe(Ts) = 0
    Substitute Xx for all Ts
    ############################################################################################## */
double BisectionZeroFunction(LocalData& dat, double Xx);
void BisectionEnergyCalc (LocalData& dat);


/*   FUNCTION TO SOLVE ENERGY BALANCE FOR MELT CONDITIONS (Qm)
     SOLVE ENERGY EQUATION
     (1-albedo)QswIn + Qlwin + QlwOut(Ts) + Qh(Ts) + Qc(Ts) +  Qe(Ts) = Qm
     Substitute Xx for all Ts
     ############################################################################################### */
void MeltEnergyCalc (LocalData& dat);


/*   FUNCTION TO SOLVE ENERGY BALANCE FOR NO SNOW CONDITIONS (Zs==0)
     SOLVE ENERGY EQUATION
     (1-albedo)QswIn + Qlwin + QlwOut(Ts) + Qh(Ts) +  Qe(Ts) = Qex
     #############   Qex is the Qc term   ##################
     ############################################################################################### */
void GroundEnergyCalc (LocalData& dat);


//  FUNCTION TO CALCULATE MELT & SUBLIMATION RATE WHEN SNOW IS PRESSENT
void MeltSublRateCalc (LocalData& dat);


//  FUNCTION TO CALCULATE MELT & SUBLIMATION RATE WHEN *NO* SNOW IS PRESSENT
void EvapCalc (EnergyBalance& eb);


/*  ####     FUNCTION TO  CALCULATE THE CHANGE IN SNOWPACK DEPTH
    All Input variables are in (SWE) Snow Water equivialnce ***Even Ps (Precipitation as Snow) In the input files ***
    Here Ps gets converted to Snowpack dpeth and stays that way
    Zs, DELTz, TotwLoss is the Snowpack depth, ***NOT IN SWE EVER***
    To Convert from SWE to Snow pack depth (Zs) Multiply by density of water (ROWw) and Devide by Density of Snow (ROWs) or (ROWFs)
    SublL, Ml, Are always in SWE
*/
void DeltaSnowPack (EnergyBalance& eb);


// FUNCTION TO MAKE SURE PROPER MASS OF WATER IS DELEVERED TO ATS WHEN SNOWPACK DISAPEARS
void WaterMassCorr(EnergyBalance& eb);


/*   FUNCTION TO CALCULATE ENERGY BALANCE AND Ml WHEN SNOWPACK IS TEENY TINY
     This section is for to trim the enery deleiverd/taken from the soil with the snowpack is really really small
     Qc = -Ks*(Ts-Tb)/Zs;  --> blows up.  Boom!
     When Zs < 0.009 [m] or 1 cm we callcualte energy balance fromt he soil surface
     *****  This is ONLY done with the snow is melting ******
     **  Otherwise Snow Energy Balance is calucalted ** */
void TeenyTinySnowPack (LocalData& tiny);


// FUNCTION TO ADDS UP ALL THE CHANGES TO THE SNOWPACK
void SnowPackCalc (EnergyBalance& eb);


// FUNCTION TO TRACKS THE TIME (IN DAYS) WHERE NO NEW SNOW AS FALLEN ~> USED IN SNOW DENSITY
void TrackSnowDays (EnergyBalance& eb);


// FUNCTION TO CALCULATE SNOWPACK DENSITY ~> WEIGHTED AVERAGE OVER THREE POTEINTAL LAYERS OF SNOW
void SnowPackDensity (EnergyBalance& eb);


// FUNCTION TO CONVERT TO SWE TO FRESHLY FALLEN SNOW DEPTH
void SWE (EnergyBalance& eb);


// FUNCTION TO CALCULATE WATER TEMPURATER
void WaterTemp (EnergyBalance& eb);

// MAIN SNOW ENERGY BALANCE FUNCTION
void SnowEnergyBalance (LocalData& seb);

}// Namespace

#endif
