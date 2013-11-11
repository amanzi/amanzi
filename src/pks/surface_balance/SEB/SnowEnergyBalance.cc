#include <iostream>
#include <cmath>

//#include "dbc.hh"
#include "SnowEnergyBalance.hh"

/* THIS CODE WILL CALCULATE THE SNOWSURFACE ENERGY BALANCE THE SNOW SURFACE TEMPERUATURE.

**** Incomming Longwave radation is cacualted in this version, but if data is avaialbe we could incoroerate it with the avialable met data****

+++Atmospheric pressure is often used in snow models, If data is available we could incorperate it but for now Pa is held constant at
100 pa++++

*** Equation for saturated vapor pressure over water is taken from Bolton, 1980 'Monthly Weather Review'
*** Equation for saturated vaport pressure over snow is taken from Buck, 1996 'Buck Research Manual'
See: http://cires.colorado.edu/~voemel/vp.html
**********************   */


void SurfaceEnergyBalance::CalcEFluxTempIndependent (LocalData& seb) {
  // Calculate incoming long-wave radiation
  seb.st_energy.fQlwIn = 1.08*(1-std::exp(-0.01*(std::pow(std::pow(10,11.4-(2353/seb.vp_air.dewpoint_temp)),(seb.st_energy.air_temp/2016)))))*seb.st_energy.stephB*(std::pow(seb.st_energy.air_temp,4));

  // Calculate D_h, D_e
  seb.st_energy.Dhe=(((std::pow(seb.st_energy.VKc,2)*seb.st_energy.Us))/(std::pow(log(seb.st_energy.Zr/seb.st_energy.Zo),2)));
}


void SurfaceEnergyBalance::CalcEFluxTempDependent (LocalData& seb, double T) {
  double Sqig;
  if (seb.st_energy.Us == 0.) {
    Sqig = 0.;
  } else {
    double Ri  = ((seb.st_energy.gZr*(seb.st_energy.air_temp-T))/(seb.st_energy.air_temp*std::pow(seb.st_energy.Us,2)));
    Sqig = (1/(1+10*Ri));
  }

  // Calculate outgoing long-wave radiation
  seb.st_energy.fQlwOut = -seb.st_energy.SEs*seb.st_energy.stephB*(std::pow(T,4));

  // Calculate sensible heat flux
  seb.st_energy.fQh = seb.st_energy.rowaCp*seb.st_energy.Dhe*Sqig*(seb.st_energy.air_temp-T);

  // Update vapore pressure of snow
  seb.vp_snow.temp=T;
  VaporCalc (seb.vp_snow);

  // Calculate latent heat flux
  seb.st_energy.fQe = seb.st_energy.rowaLs*seb.st_energy.Dhe*Sqig*0.622*((seb.vp_air.actual_vaporpressure-seb.vp_snow.saturated_vaporpressure)/seb.st_energy.Apa);

  // Calculate heat conducted to ground
       // if (seb.st_energy.ht_snow<=0.009) { // if Zs really really small Qc = -Ks*(Ts-Tb)/Zs;  --> blows up
       //     seb.st_energy.fQc = seb.st_energy.snow_conduct_value*(T-seb.st_energy.Tb)/seb.st_energy.snow_groundTrans;
       // }else{
       //seb.st_energy.fQc = seb.st_energy.snow_conduct_value*(T-seb.st_energy.Tb)/seb.st_energy.ht_snow;
       //  }
    seb.st_energy.fQc = seb.st_energy.snow_conduct_value*(T-seb.st_energy.Tb);
}


// #### FUNCTIONS TO CALCULATE SATURATION VAPOR PRESSURE & ACTUALL VAPOR PRESSURE #########################
void SurfaceEnergyBalance::VaporCalc (VaporPressure& vp) {
  double temp;
  //Convert from Kelvin to Celsius
  temp = vp.temp-273.15;
  // Sat vap. press o/water Dingman D-7 (Bolton, 1980)
  vp.saturated_vaporpressure = 0.611*std::exp((17.67*temp)/(temp+243.5));
  // (Bolton, 1980)
  vp.actual_vaporpressure=vp.saturated_vaporpressure*vp.relative_humidity;
  // Find dewpoint Temp Dingman D-11
  vp.dewpoint_temp=(std::log(vp.actual_vaporpressure)+0.4926)/(0.0708-0.00421*std::log(vp.actual_vaporpressure));
  // Convert Tdp from Celsius to Kelvin
  vp.dewpoint_temp=vp.dewpoint_temp + 273.15;
}


// #### FUNCTIONS TO CALCULATE ALBEDO
//void SurfaceEnergyBalance::AlbedoCalc (EnergyBalance& eb) {
//  // If snow is present
//  if (eb.ht_snow>0.0) {
//    // Albedo -- Ling and Zhang (2004) Method
//    if (eb.density_snow<=432.238) {
//      eb.albedo_value=1.0-0.247*(std::pow((0.16+110*std::pow((eb.density_snow/1000),4)),0.5));
//    } else {
//      eb.albedo_value=0.6-(eb.density_snow/4600);
//    }
//  } else {
//    //If no snow is present Calculate Surface Albedo
//    if (eb.water_depth>0.0) {// Checking for standing water oughta be deeper then the plants that stand say 10 cm
//      eb.albedo_value=0.6; // ******* Refine the albedo calculation using clm *****
//    }
//    if (eb.water_depth<=0.0) {
//      eb.albedo_value=0.15; // Albedo for Tundra heather  Dingman Table D-2
//    }
//  }
//}

// #### FUNCTIONS TO CALCULATE ALBEDO
void SurfaceEnergyBalance::AlbedoCalc (EnergyBalance& eb) {
    double perSnow = 0.0, perTundra=0.0, perWater=0.0;
    double AlTundra=0.15, AlWater=0.6;
    double AlSnow = 0.0;
    double TransitionVal=eb.AlbedoTrans;  // Set to 2 cm

    if (eb.density_snow<=432.238) {
        AlSnow=1.0-0.247*(std::pow((0.16+110*std::pow((eb.density_snow/1000),4)),0.5));
    } else {
        AlSnow=0.6-(eb.density_snow/4600);
    }
    if (eb.ht_snow>TransitionVal) {
        perSnow = 1;  // Snow is too deep for albedo wt ave
    }else{
        if (eb.water_depth<=0.0) {  // dry ground
            perSnow = pow((eb.ht_snow/TransitionVal),2); // Transition to dry ground
            perTundra = 1-perSnow;
        }else {
            if (eb.ht_snow>TransitionVal) { // Transition for surface water albedo = 1cm
                perSnow = 1;  // Snow is too deep over surface water for albedo wt ave
            }else{
               perSnow = eb.ht_snow/TransitionVal;  // Transitions to surface water
               perWater = 1-perSnow;
            }
        }
    }
    // weighted vaverage function for surface albedo
    eb.albedo_value = ((AlSnow*perSnow)+(AlTundra*perTundra)+(AlWater*perWater))/(perSnow + perTundra + perWater);
}



// ### FUNCTION TO CALCULATE THERMAL CONDUCTIVITY OF SNOW
void SurfaceEnergyBalance::ThermalConductSnow (EnergyBalance& eb) {
  // If snow is present
  if (eb.ht_snow>0.0) {
    eb.snow_conduct_value=2.9e-6*std::pow(eb.density_snow,2);
  } else {
    eb.snow_conduct_value=0;
  }
}


/*  FUNCTION TO SOLVE SNOW SURFACE ENERGY BALANCE FOR SURFACE TEMP (Ts)
    ###############################################################################################*
    Bisection Method
    (1-albedo)QswIn + Qlwin + QlwOut(Ts) + Qh(Ts) + Qc(Ts) +  Qe(Ts) = 0
    Substitute Xx for all Ts
    ####################################################################################################
*/
double SurfaceEnergyBalance::BisectionZeroFunction(LocalData& seb, double Xx) {
  CalcEFluxTempDependent(seb, Xx);

  // BALANCE EQUATION
  double ZERO = seb.st_energy.ht_snow * (seb.st_energy.fQswIn + seb.st_energy.fQlwIn + seb.st_energy.fQlwOut + seb.st_energy.fQh
      + seb.st_energy.fQe) - seb.st_energy.fQc;
  return ZERO;
}


void SurfaceEnergyBalance::BisectionEnergyCalc (LocalData& seb) {
  double tol = 1.e-6;
  double deltaX = 5;

  double Xx = seb.st_energy.air_temp;
  double FXx = BisectionZeroFunction(seb, Xx);
  // NOTE: decreasing function
  double a,b,Fa,Fb;
  if (FXx > 0) {
    b = Xx;
    Fb = FXx;
    a = Xx;
    Fa = FXx;
    while (Fa > 0) {
      b = a;
      Fb = Fa;
      a += deltaX;
      Fa = BisectionZeroFunction(seb,a);
    }
  } else {
    a = Xx;
    Fa = FXx;
    b = Xx;
    Fb = FXx;
    while (Fb < 0) {
      a = b;
      Fa = Fb;
      b -= deltaX;
      Fb = BisectionZeroFunction(seb,b);
    }
  }

 // ASSERT(Fa*Fb < 0);

  int maxIterations = 200;
  double res;
  int iter;
  for (int i=0; i<maxIterations; ++i) { //Besection Iterations Loop Solve for Ts using Energy balance equation #######
    Xx = (a+b)/2;
    res = BisectionZeroFunction(seb, Xx);

    if (res>0) {
      b=Xx;
    } else {
      a=Xx;
    }
    if (std::abs(res)<tol) {
      break;
    }
    iter=i;
  }//End Bisection Interatiion Loop  Solve for Ts using Energy balance equation ##################################

  if (std::abs(res) >= tol) {
  //  ASSERT(0);
  }
  seb.st_energy.Ts=Xx;
}


/*   FUNCTION TO SOLVE ENERGY BALANCE FOR MELT CONDITIONS (Qm)
     SOLVE ENERGY EQUATION
     (1-albedo)QswIn + Qlwin + QlwOut(Ts) + Qh(Ts) + Qc(Ts) +  Qe(Ts) = Qm
     Substitute Xx for all Ts
     ####################################################################################################
*/
void SurfaceEnergyBalance::MeltEnergyCalc (LocalData& seb) {
  CalcEFluxTempDependent(seb, seb.st_energy.Ts);

  // Melt energy is the balance
  seb.st_energy.Qm = seb.st_energy.ht_snow * (seb.st_energy.fQswIn + seb.st_energy.fQlwIn + seb.st_energy.fQlwOut
      + seb.st_energy.fQh + seb.st_energy.fQe) - seb.st_energy.fQc;
    
    if (seb.st_energy.ht_snow<=seb.st_energy.snow_groundTrans) { // if Zs really really small Qc = -Ks*(Ts-Tb)/Zs;  --> blows up
         seb.st_energy.Qm = seb.st_energy.Qm/seb.st_energy.snow_groundTrans;
     }else{
    seb.st_energy.Qm = seb.st_energy.Qm/seb.st_energy.ht_snow;
      }
}


/*   FUNCTION TO SOLVE ENERGY BALANCE FOR NO SNOW CONDITIONS (Zs==0)
     SOLVE ENERGY EQUATION
     (1-albedo)QswIn + Qlwin + QlwOut(Ts) + Qh(Ts) +  Qe(Ts) = Qc
     #############   Qex is the Qc term   ##################
     ####################################################################################################
*/
void SurfaceEnergyBalance::GroundEnergyCalc (LocalData& seb) {
  seb.st_energy.fQlwOut = -seb.st_energy.SEs*seb.st_energy.stephB*(std::pow(seb.st_energy.Tb,4));

  double Sqig;
  if (seb.st_energy.Us == 0.) {
    Sqig = 0.;
  } else {
    double Ri = ((seb.st_energy.gZr*(seb.st_energy.air_temp-seb.st_energy.Tb))/(seb.st_energy.air_temp*std::pow(seb.st_energy.Us,2)));
      if (Ri<0) {// Unstable condition
          Sqig = (1-10*Ri);
      }else{// Stable Condition
    Sqig = (1/(1+10*Ri));
      }
  }

  seb.st_energy.fQh = seb.st_energy.rowaCp*seb.st_energy.Dhe*Sqig*(seb.st_energy.air_temp-seb.st_energy.Tb);

  if (seb.st_energy.water_depth>0.0) {// Checking for standing water
    VaporCalc (seb.vp_ground);
    seb.st_energy.fQe = seb.st_energy.rowaLe*seb.st_energy.Dhe*Sqig*0.622*((seb.vp_air.actual_vaporpressure-seb.vp_ground.saturated_vaporpressure)/seb.st_energy.Apa);
  } else {// no standing water
    seb.st_energy.fQe = seb.st_energy.porrowaLe*seb.st_energy.Dhe*Sqig*0.622*((seb.vp_air.actual_vaporpressure-seb.vp_ground.actual_vaporpressure)/seb.st_energy.Apa);
  }

  // Heat flux to ground surface is the balance.
  seb.st_energy.fQc=seb.st_energy.fQswIn+seb.st_energy.fQlwIn+seb.st_energy.fQlwOut+seb.st_energy.fQh+seb.st_energy.fQe;
}


// FUNCTION TO CALCULATE Qc CONDUCTIVE HEAT FLUX THROUGH THE SNOW PACK
// This value needs to be calculated here because Zs was multiblied through Energy balance equation
// in order to avoid    Qc = -Ks*(Ts-Tb)/Zs;  --> blowing up when Zs is really small
void SurfaceEnergyBalance::CalcQc (LocalData& seb){
    // Calculate heat conducted to ground
     if (seb.st_energy.ht_snow<=seb.st_energy.snow_groundTrans) { // if Zs really really small Qc = -Ks*(Ts-Tb)/Zs;  --> blows up
         seb.st_energy.fQc = seb.st_energy.snow_conduct_value*(seb.st_energy.Ts-seb.st_energy.Tb)/seb.st_energy.snow_groundTrans;
     }else{
    seb.st_energy.fQc = seb.st_energy.snow_conduct_value*(seb.st_energy.Ts-seb.st_energy.Tb)/seb.st_energy.ht_snow;
      }
}


//  FUNCTION TO CALCULATE MELT & SUBLIMATION RATE WHEN SNOW IS PRESSENT
void SurfaceEnergyBalance::MeltSublRateCalc (LocalData& seb) {
  // Calculate water melted  *** Equation from UEB (49)
  seb.st_energy.Mr=seb.st_energy.Qm/(seb.st_energy.density_w*seb.st_energy.Hf);      // Change Mr = Qm/(ROWw*Hf) to --> Mr = Qm/(ROWw*Hf) + (Pr/Dt); ** this will mean changing DeltaSnowPack & WaterMassCorr
  seb.st_energy.Ml=seb.st_energy.Mr*seb.st_energy.Dt;                                        // Ml=Mr*Dt;
  seb.st_energy.SublR=-seb.st_energy.fQe/(seb.st_energy.density_w*seb.st_energy.Ls); // SublR=-fQe/(ROWw*Ls);  // SublR is a rate [m/s]
  seb.st_energy.SublL=seb.st_energy.SublR*seb.st_energy.Dt;                                  // SublL=SublR*Dt; // SublL is a swe lenght [m]
}


//  FUNCTION TO CALCULATE MELT & SUBLIMATION RATE WHEN *NO* SNOW IS PRESSENT
void SurfaceEnergyBalance::EvapCalc (EnergyBalance& eb) {
  eb.EvR=eb.fQe/(eb.density_w*eb.Le); // EvR=Qe/(ROWw*Le); // SublR is a rate [m/s]
  eb.EvL=eb.EvR*eb.Dt;      // EvL=EvR*Dt;           // SublL is a swe lenght [m]
  if (eb.ht_snow == 0.) {// Makesure Ml is not updated when a tiny amount of snow is present, Already taken care of by MeltSublRateCalc
    eb.Ml = eb.Pr+eb.EvL;  // Ml = Pr+EvL;
    eb.Mr=eb.Ml/eb.Dt;     // Mr=Ml/Dt;
  }
}


/*  ####     FUNCTION TO  CALCULATE THE CHANGE IN SNOWPACK DEPTH
    All Input variables are in (SWE) Snow Water equivialnce ***Even Ps (Precipitation as Snow) In the input files ***
    Here Ps gets converted to Snowpack dpeth and stays that way
    Zs, DELTz, TotwLoss is the Snowpack depth, ***NOT IN SWE EVER***
    To Convert from SWE to Snow pack depth (Zs) Multiply by density of water (ROWw) and Devide by Density of Snow (ROWs) or (ROWFs)
    SublL, Ml, Are always in SWE
*/
void SurfaceEnergyBalance::DeltaSnowPack (EnergyBalance& eb) {
  //if (SublL<0) {    //briefly changing SublL into snow (frost) lenght units *** This is for condensation ****
  //    CiL=-SublL*ROWw/ROWfrost, SublL=0 };  // I assume that ice is condensing on the snow surface therefore density is near that of ice
  //if (SublL>=0) {   //briefly changing SublL into snow lengh units *** This is for sublimation ****
  //    SublL=SublL*ROWw/ROWs };
  //TotwLoss = Ml*(ROWw/ROWs) + SublL - CiL;  //Calculate change in snow pack
  //if (SublL>=0) {  SublL=SublL*ROWs/ROWw };                     //Changing SublL back to SWE *** This is for sublimation ****
  if (eb.SublL<0) {    //briefly changing SublL into snow (frost) lenght units *** This is for condensation ****
    eb.CiL = -eb.SublL * eb.density_w / eb.density_frost;  // Assume ice is condensing on snow surface with a density near ice ~> 800[kg/m^3]
    eb.SublL = 0;
  } else {   //briefly changing SublL into snow lengh units *** This is for sublimation ****
    eb.SublL = eb.SublL * eb.density_w / eb.density_snow;
    eb.CiL = 0.;
  }

  eb.TotwLoss = eb.Ml * (eb.density_w / eb.density_snow) + eb.SublL - eb.CiL;  //Calculate change in snow pack

  if (eb.SublL >= 0) { //Changing SublL back to SWE *** This is for sublimation ****
    eb.SublL = eb.SublL * eb.density_snow / eb.density_w;
  }
}


// FUNCTION TO MAKE SURE PROPER MASS OF WATER IS DELEVERED TO ATS WHEN SNOWPACK DISAPEARS
void SurfaceEnergyBalance::WaterMassCorr(EnergyBalance& eb) {
  //if (Ml>0) {
  //   if (SublL>0) {                             //When SublL is great then 0 Sublimation! Find the ratio of water atributed to SublL and subtract it from the water delivered to ATS
  //      SublL=(SublL*(ROWw/ROWs)/abs(TotwLoss));//SublL is temporarily changed to snow length
  //      SublL=SublL*Zs*(ROWs/ROWw);             //SublL is changed back to SWE & ratio of avaialbe snowpack for sublimation is found
  //      Ml=(Zs+Ps+(SublL*(ROWw/ROWs)))*(ROWs/ROWw);
  //   }else{Ml=(Zs+Ps+(SublL*-1))*(ROWs/ROWw) }; //SublL is less then 0, Frost! Add it to the water delivered to ATS
  //if (Ml==0) { SublL=(Zs+Ps)*(ROWs/ROWw) };     // The snowpack sublimated away.
  //TotwLoss=Zs;
  //Mr=Ml/Dt;}
  double Tswe = eb.density_w/eb.density_snow;
  double TZs = eb.density_snow/eb.density_w;
  if (eb.Ml>0) {
    if (eb.SublL>0) {               // When SublL is great then 0 Sublimation! Find the ratio of water atributed to SublL and subtract it from the water delivered to ATS
      eb.SublL=(eb.SublL*(Tswe)/std::abs(eb.TotwLoss));  //SublL is temporarily changed to snow length
      eb.SublL=eb.SublL*eb.ht_snow*(TZs);  //SublL is changed back to SWE & ratio of avaialbe snowpack for sublimation is found
      eb.Ml=(eb.ht_snow+eb.Ps+(eb.SublL*(Tswe)))*(TZs);
    } else {
      eb.Ml=(eb.ht_snow+eb.Ps+(eb.SublL*-1))*(TZs); //SublL is less then 0, Frost! Add it to water delivered to ATS
    }
  }
  if (eb.Ml==0) {
    eb.SublL=(eb.ht_snow+eb.Ps)*(TZs); // The snowpack sublimated away.
  }
  eb.TotwLoss = eb.ht_snow;
  eb.Mr = eb.Ml/eb.Dt;
}



/*   FUNCTION TO CALCULATE ENERGY BALANCE AND Ml WHEN SNOWPACK IS TEENY TINY
     This section is for to trim the enery deleiverd/taken from the soil with the snowpack is really really small
     Qc = -Ks*(Ts-Tb)/Zs;  --> blows up.  Boom!
     When Zs < 0.009 [m] or 1 cm we callcualte energy balance fromt he soil surface
     *****  This is ONLY done with the snow is melting ******
     **  Otherwise Snow Energy Balance is calucalted **
        ***** THIS IS NOW TURNED OFF ****** 
 BECAUSE IT CAN NOT HANDLE WHEN SNOW ACCUMULATION IS REALLY SMALL
 WHEN FOR EXAMPLE THE TIME STEP BECOMES REALLY SMALL ON A SNOW DAY
   Maybe fix this if you got nothing better to do.
 Now Zs can be no smaller then 0.009 *just for caluclating Qc*
     */
//void SurfaceEnergyBalance::TeenyTinySnowPack (LocalData& seb) {
//  seb.st_energy.Qm=seb.st_energy.Mr*(seb.st_energy.density_w*seb.st_energy.Hf); //Recalculating Qm based off of Snowpack limit
//  double ZsHold = seb.st_energy.ht_snow;
//  seb.st_energy.ht_snow=0.0;
//  AlbedoCalc(seb.st_energy); // Calculating Surface Albedo (No Snow)
//  seb.st_energy.fQswIn=(1-seb.st_energy.albedo_value)*seb.st_energy.QswIn;
//  GroundEnergyCalc(seb); //  Solving Energy balance for the ground
//  seb.st_energy.funcall="BARE-teenyseb";
//  seb.st_energy.ht_snow=ZsHold;
//  EvapCalc(seb.st_energy);//Calculating Evaporation from bare-ground "ht_snow" is keyed to *NOT* recalculate Ml because WaterMassCorr has already done that
//  seb.st_energy.Ml=seb.st_energy.Ml+seb.st_energy.Pr+seb.st_energy.EvL;
//}


// FUNCTION TO ADDS UP ALL THE CHANGES TO THE SNOWPACK
void SurfaceEnergyBalance::SnowPackCalc (EnergyBalance& eb) {
  // Zs = Zs + Ps - TotwLoss;    // New formate
    
    if (eb.ht_snow>0) {
  eb.ht_snow = eb.ht_Zs_settled + eb.Ps - eb.TotwLoss;  //DELTz is the new snowpack depth
    }else{
      eb.ht_snow = eb.ht_snow + eb.Ps - eb.TotwLoss;  
    }
//  ASSERT(eb.ht_snow >= 0.);
}


// FUNCTION TO TRACKS THE TIME (IN DAYS) WHERE NO NEW SNOW AS FALLEN ~> USED IN SNOW DENSITY
void SurfaceEnergyBalance::TrackSnowDays (EnergyBalance& eb) {
    double Beta=((eb.ht_snow*eb.density_snow)+(eb.Ps*eb.density_freshsnow)+(eb.CiL*eb.density_frost));
    if (eb.ht_snow>0) {        
        eb.nosnowdays=(((eb.nosnowdays + (eb.Dt / 86400))*eb.ht_snow*eb.density_snow) + ((eb.Dt / 86400)*eb.Ps*eb.density_freshsnow)+((eb.NDSfrost + (eb.Dt / 86400))*eb.CiL*eb.density_frost))/Beta;
    }else{
     eb.nosnowdays = 0;   
    }
  //if (eb.Ps < 0.0001) {//If less then a mm of snow
  //  eb.nosnowdays += eb.Dt / 86400;
  //} else {
  //  eb.nosnowdays = 0;
 // }
  }

// CALCULATES SNOW DEFORMATION ~> NEW DENSITY AND HEIGHT OFF AGED SNOW 'LAYER' (Martinec, 1977)
void SurfaceEnergyBalance::SnowDeformationModel (EnergyBalance& eb) {
    // Track days with now snow and formulate snow deformation
    double ndensity = std::pow((eb.nosnowdays+1),0.3);
    eb.density_snow=eb.density_freshsnow*ndensity;// New Density of settled snow
    eb.ht_Zs_settled = eb.ht_snow * (eb.HoldROWs / eb.density_snow); // New Height of settled snow
}

// FUNCTION TO CALCULATE SNOWPACK DENSITY ~> WEIGHTED AVERAGE OVER THREE POTEINTAL LAYERS OF SNOW
void SurfaceEnergyBalance::SnowPackDensity (EnergyBalance& eb) {
  if (eb.ht_snow == 0.) {
    eb.density_snow=100;
    eb.HoldROWs=100;
  } else {
    //Weighted average of the three layers of the snowpack
    //ROWs=((Ps*ROWfs)+(Zs_settled*ROWs)+(CiL*ROWfrost))/(Ps+Zs_settled+CiL);  // Wt. average for snow density
    double denominator=eb.Ps+eb.ht_Zs_settled+eb.CiL;
    eb.density_snow=((eb.Ps*eb.density_freshsnow)+(eb.ht_Zs_settled*eb.density_snow)+(eb.CiL*eb.density_frost))/(denominator);
    if (eb.density_snow>950) {// Capping snow density ~> it should NEVER get this dense anyway.
      eb.density_snow=950.1;
    }
  }
}


// FUNCTION TO CONVERT TO SWE TO FRESHLY FALLEN SNOW DEPTH
void SurfaceEnergyBalance::SWE (EnergyBalance& eb) {
  eb.Ps = eb.Ps*(eb.density_w/eb.density_freshsnow);
}

// FUNCTION TO FIND TEMPURATURE OF WATER FLOWING INTP ATS
void SurfaceEnergyBalance::WaterTemp(EnergyBalance& eb) {
  if (eb.Ml<=0.0) {
    eb.Trw=eb.air_temp;
  } else {
    eb.Trw=273.15;
  }
}


// MAIN SNOW ENERGY BALANCE FUNCTION
void SurfaceEnergyBalance::SnowEnergyBalance (LocalData& seb) {
  // Caculate Vapor pressure and dewpoint temperature from Air
  VaporCalc(seb.vp_air);

  // Find Albedo
  AlbedoCalc(seb.st_energy); // Calculating Surface Albedo
  seb.st_energy.fQswIn=(1-seb.st_energy.albedo_value)*seb.st_energy.QswIn;

  // Find Thermal Conductivity of snow
  ThermalConductSnow(seb.st_energy); //Calculating themeral conductifity of snow

  // Convert snow leath Ps (m) in SWE lenght
  SWE(seb.st_energy);

  // Update temperature-independent fluxes
  CalcEFluxTempIndependent(seb);

  if (seb.st_energy.ht_snow>0) { // If snow
    //Solving for snowsurface temperature and Energy balance equation
    BisectionEnergyCalc(seb);
    seb.st_energy.funcall="SNOW";
    if (seb.st_energy.Ts<=273.15) { // Snow is not melting
      seb.st_energy.Qm=0; //  no water leaving snowpack as melt water
      seb.st_energy.Ml=0;
    } else {
      seb.st_energy.Ts=273.15; // Setting snow temperature to zero
      MeltEnergyCalc(seb); // Recaculating Energy Balance when melt in occuring
      seb.st_energy.funcall="MELT";
    }

    // Calculating melt, sublimation, and condensation rate.
    MeltSublRateCalc(seb);

    // Calculate the change in snowpack depth.
    DeltaSnowPack(seb.st_energy);

    // Make sure proper mass of snowpack water gets delivered to AT
    if (seb.st_energy.ht_snow < ((seb.st_energy.Ps - seb.st_energy.TotwLoss)*-1)) {
      WaterMassCorr(seb.st_energy);
      seb.st_energy.funcall="Melt too small";
    }
      
    // Rain on Snow = Rain water flows through snowpack added to water delivered to ATS & Assinged temperarture of 0 Celcius
      seb.st_energy.Mr = seb.st_energy.Mr + (seb.st_energy.Pr / seb.st_energy.Dt);
      
    // Correct Qc & check for small snowpack
      CalcQc(seb);
      
    // Check for small melting snowpack, which can cause super high Qc:
    // Qc = -Ks*(Ts-Tb)/Zs;  --> blows up
    //if (seb.st_energy.ht_snow<=0.009) {
    //  if (seb.st_energy.Ts==273.15) {
    //    TeenyTinySnowPack(seb);
    //  }
    //}

  } else { // no snow
    GroundEnergyCalc(seb);
    seb.st_energy.funcall="BARE";
    EvapCalc(seb.st_energy);
    seb.st_energy.TotwLoss = 0.; // no snow to melt
  }
    
  //Calculate Snow Deformation
  seb.st_energy.HoldROWs=seb.st_energy.density_snow;
  SnowDeformationModel(seb.st_energy);

  // Update snowpack density
  SnowPackDensity(seb.st_energy);

  // Update snow height
  SnowPackCalc(seb.st_energy);
  if (seb.st_energy.ht_snow<=0.0000001) {
    seb.st_energy.ht_snow=0.0;
  }
  
  // Update days since last snow for snow pack deformation
  TrackSnowDays(seb.st_energy);

  // Calculate water temperature
  WaterTemp(seb.st_energy);

}
