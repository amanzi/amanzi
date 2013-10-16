#include <iostream>
#include <cmath>

#include "SnowEnergyBalance.hh"

/* THIS CODE WILL CALCULATE THE SNOWSURFACE ENERGY BALANCE THE SNOW SURFACE TEMPERUATURE.

**** Incomming Longwave radation is cacualted in this version, but if data is avaialbe we could incoroerate it with the avialable met data****

+++Atmospheric pressure is often used in snow models, If data is available we could incorperate it but for now Pa is held constant at
100 pa++++

*** Equation for saturated vapor pressure over water is taken from Bolton, 1980 'Monthly Weather Review'
*** Equation for saturated vaport pressure over snow is taken from Buck, 1996 'Buck Research Manual'
See: http://cires.colorado.edu/~voemel/vp.html
**********************   */





// #### FUNCTIONS TO CALCULATE SATURATION VAPOR PRESSURE & ACTUALL VAPOR PRESSURE #########################
void SurfaceEnergyBalance::VaporCalc (VaporPressure& vaporpressure) {
  double temp;
  //Convert from Kelvin to Celsius
  temp = vaporpressure.temp-273.15;
  // Sat vap. press o/water Dingman D-7 (Bolton, 1980)
  vaporpressure.saturated_vaporpressure = 0.611*std::exp((17.67*temp)/(temp+243.5));
  // (Bolton, 1980)
  vaporpressure.actual_vaporpressure=vaporpressure.saturated_vaporpressure*vaporpressure.relative_humidity;
  // Find dewpoint Temp Dingman D-11
  vaporpressure.dewpoint_temp=(std::log(vaporpressure.actual_vaporpressure)+0.4926)/(0.0708-0.00421*std::log(vaporpressure.actual_vaporpressure));
  // Convert Tdp from Celsius to Kelvin
  vaporpressure.dewpoint_temp=vaporpressure.dewpoint_temp + 273.15;
}


// #### FUNCTIONS TO CALCULATE ALBEDO
void SurfaceEnergyBalance::AlbedoCalc (LocalData& albedo) {
  // If snow is present
  if (albedo.st_energy.ht_snow>0.0) {
    // Albedo -- Ling and Zhang (2004) Method
    if (albedo.st_energy.density_snow<=450) {
      albedo.st_energy.albedo_value=1.0-0.247*(std::pow((0.16+110*std::pow((albedo.st_energy.density_snow/1000),4)),0.5));
    } else {
      albedo.st_energy.albedo_value=0.6-(albedo.st_energy.density_snow/4600);
    }
  } else {
    //If no snow is present Calculate Surface Albedo
    if (albedo.st_energy.water_depth>0.0) {// Checking for standing water oughta be deeper then the plants that stand say 10 cm
      albedo.st_energy.albedo_value=0.6; // ******* Refine the albedo calculation using clm *****
    }
    if (albedo.st_energy.water_depth<=0.0) {
      albedo.st_energy.albedo_value=0.15; // Albedo for Tundra heather  Dingman Table D-2
    }
  }
}



// ### FUNCTION TO CALCULATE THERMAL CONDUCTIVITY OF SNOW
void SurfaceEnergyBalance::ThermalConductSnow (LocalData& thermalsnow) {
  // If snow is present
  if (thermalsnow.st_energy.ht_snow>0.0) {
    thermalsnow.st_energy.snow_conduct_value=2.9e-6*std::pow(thermalsnow.st_energy.density_snow,2);
  } else {
    thermalsnow.st_energy.snow_conduct_value=0;
  }
}


/*  FUNCTION TO SOLVE SNOW SURFACE ENERGY BALANCE FOR SURFACE TEMP (Ts)
    ###############################################################################################*
    Bisection Method
    (1-albedo)QswIn + Qlwin + QlwOut(Ts) + Qh(Ts) + Qc(Ts) +  Qe(Ts) = 0
    Substitute Xx for all Ts
    ####################################################################################################
*/
void SurfaceEnergyBalance::BisectionEnergyCalc (LocalData& snow) {
  double a = snow.st_energy.air_temp+20; //Setting bounds of bisection proximal to Air Temp
  double b = snow.st_energy.air_temp-20; //Setting bounds of bisection proximal to Air Temp
  double Xx =0.0;
  double Ri=0.0;
  double Sqig=0.0;
  double ZERO=0.0;
  double tol=1e-8;
  int Iterations=100;
  int iter=0;

  snow.st_energy.Dhe=(((std::pow(snow.st_energy.VKc,2)*snow.st_energy.Us))/(std::pow(log(snow.st_energy.Zr/snow.st_energy.Zo),2)));
  snow.st_energy.fQlwIn = 1.08*(1-std::exp(-0.01*(std::pow(std::pow(10,11.4-(2353/snow.vp_air.dewpoint_temp)),(snow.st_energy.air_temp/2016)))))*snow.st_energy.stephB*(std::pow(snow.st_energy.air_temp,4));

  for (int i=0; i<Iterations; i++) { //Besection Iterations Loop Solve for Ts using Energy balance equation #######
    Xx = (a+b)/2;
    snow.st_energy.fQlwOut = -snow.st_energy.SEs*snow.st_energy.stephB*(std::pow(Xx,4));
    Ri = ((snow.st_energy.gZr*(snow.st_energy.air_temp-Xx))/(snow.st_energy.air_temp*std::pow(snow.st_energy.Us,2)));
    Sqig = (1/(1+10*Ri));
    snow.st_energy.fQh = snow.st_energy.rowaCp*snow.st_energy.Dhe*Sqig*(snow.st_energy.air_temp-Xx);
    snow.vp_snow.temp=Xx;
    VaporCalc (snow.vp_snow);
    snow.st_energy.fQe = snow.st_energy.rowaLs*snow.st_energy.Dhe*Sqig*0.622*((snow.vp_air.actual_vaporpressure-snow.vp_snow.saturated_vaporpressure)/snow.st_energy.Apa);
    snow.st_energy.fQc = snow.st_energy.snow_conduct_value*(Xx-snow.st_energy.Tb)/snow.st_energy.ht_snow;
    // BALANCE EQUATION
    ZERO = snow.st_energy.fQswIn + snow.st_energy.fQlwIn + snow.st_energy.fQlwOut + snow.st_energy.fQh + snow.st_energy.fQe - snow.st_energy.fQc;
    if (ZERO>0) {
      b=Xx;
    } else {
      a=Xx;
    }
    if (std::abs(ZERO)<tol) {
      break;
    }
    iter=i;
  }//End Bisection Interatiion Loop  Solve for Ts using Energy balance equation ##################################
  snow.st_energy.Ts=Xx;
}


/*   FUNCTION TO SOLVE ENERGY BALANCE FOR MELT CONDITIONS (Qm)
     SOLVE ENERGY EQUATION
     (1-albedo)QswIn + Qlwin + QlwOut(Ts) + Qh(Ts) + Qc(Ts) +  Qe(Ts) = Qm
     Substitute Xx for all Ts
     ####################################################################################################
*/
void SurfaceEnergyBalance::MeltEnergyCalc (LocalData& melt) {
  double Ri=0.0;
  double Sqig=0.0;
  melt.st_energy.fQlwIn = 1.08*(1-std::exp(-0.01*(std::pow(std::pow(10,11.4-(2353/melt.vp_air.dewpoint_temp)),(melt.st_energy.air_temp/2016)))))*melt.st_energy.stephB*(std::pow(melt.st_energy.air_temp,4));
  melt.st_energy.fQlwOut = -melt.st_energy.SEs*melt.st_energy.stephB*(std::pow(melt.st_energy.Ts,4));
  melt.st_energy.Dhe=(((std::pow(melt.st_energy.VKc,2)*melt.st_energy.Us))/(std::pow(log(melt.st_energy.Zr/melt.st_energy.Zo),2)));
  Ri  = ((melt.st_energy.gZr*(melt.st_energy.air_temp-melt.st_energy.Ts))/(melt.st_energy.air_temp*std::pow(melt.st_energy.Us,2)));
  Sqig = (1/(1+10*Ri));
  melt.st_energy.fQh = melt.st_energy.rowaCp*melt.st_energy.Dhe*Sqig*(melt.st_energy.air_temp-melt.st_energy.Ts);
  melt.vp_snow.temp=melt.st_energy.Ts;
  VaporCalc (melt.vp_snow);
  melt.st_energy.fQe = melt.st_energy.rowaLs*melt.st_energy.Dhe*Sqig*0.622*((melt.vp_air.actual_vaporpressure-melt.vp_snow.saturated_vaporpressure)/melt.st_energy.Apa);
  melt.st_energy.fQc = melt.st_energy.snow_conduct_value*(melt.st_energy.Ts-melt.st_energy.Tb)/melt.st_energy.ht_snow;
  melt.st_energy.Qm = melt.st_energy.fQswIn + melt.st_energy.fQlwIn + melt.st_energy.fQlwOut + melt.st_energy.fQh + melt.st_energy.fQe - melt.st_energy.fQc;

}


/*   FUNCTION TO SOLVE ENERGY BALANCE FOR NO SNOW CONDITIONS (Zs==0)
     SOLVE ENERGY EQUATION
     (1-albedo)QswIn + Qlwin + QlwOut(Ts) + Qh(Ts) +  Qe(Ts) = Qex
     #############   Qex is the Qc term   ##################
     ####################################################################################################
*/
void SurfaceEnergyBalance::GroundEnergyCalc (LocalData& grd) {
  double Ri=0.0;
  double Sqig=0.0;

  grd.st_energy.fQlwIn = 1.08*(1-std::exp(-0.01*(std::pow(std::pow(10,11.4-(2353/grd.vp_air.dewpoint_temp)),(grd.st_energy.air_temp/2016)))))*grd.st_energy.stephB*(std::pow(grd.st_energy.air_temp,4));
  grd.st_energy.fQlwOut = -grd.st_energy.SEs*grd.st_energy.stephB*(std::pow(grd.st_energy.Tb,4));
  grd.st_energy.Dhe=(((std::pow(grd.st_energy.VKc,2)*grd.st_energy.Us))/(std::pow(log(grd.st_energy.Zr/grd.st_energy.Zo),2)));
  Ri = ((grd.st_energy.gZr*(grd.st_energy.air_temp-grd.st_energy.Tb))/(grd.st_energy.air_temp*std::pow(grd.st_energy.Us,2)));
  Sqig = (1/(1+10*Ri));
  grd.st_energy.fQh = grd.st_energy.rowaCp*grd.st_energy.Dhe*Sqig*(grd.st_energy.air_temp-grd.st_energy.Tb);
  if (grd.st_energy.water_depth>0.0) {// Checking for standing water
    VaporCalc (grd.vp_ground);
    grd.st_energy.fQe = grd.st_energy.rowaLe*grd.st_energy.Dhe*Sqig*0.622*((grd.vp_air.actual_vaporpressure-grd.vp_ground.saturated_vaporpressure)/grd.st_energy.Apa);
  }
  if (grd.st_energy.water_depth<=0.0) {// ****  We thing there needs to be a porosity term in this equation MAY NEED REFINEMENT *****
    grd.st_energy.fQe = grd.st_energy.porrowaLe*grd.st_energy.Dhe*Sqig*0.622*((grd.vp_air.actual_vaporpressure-grd.vp_ground.actual_vaporpressure)/grd.st_energy.Apa);
  }
  grd.st_energy.fQc=grd.st_energy.fQswIn+grd.st_energy.fQlwIn+grd.st_energy.fQlwOut+grd.st_energy.fQh+grd.st_energy.fQe;
}


//  FUNCTION TO CALCULATE MELT & SUBLIMATION RATE WHEN SNOW IS PRESSENT
void SurfaceEnergyBalance::MeltSublRateCalc (LocalData& meltsub) {
  // Calculate water melted  *** Equation from UEB (49)
  meltsub.st_energy.Mr=meltsub.st_energy.Qm/(meltsub.st_energy.density_w*meltsub.st_energy.Hf);      // Mr=Qm/(ROWw*Hf);
  meltsub.st_energy.Ml=meltsub.st_energy.Mr*meltsub.st_energy.Dt;                                        // Ml=Mr*Dt;
  meltsub.st_energy.SublR=-meltsub.st_energy.fQe/(meltsub.st_energy.density_w*meltsub.st_energy.Ls); // SublR=-fQe/(ROWw*Ls);  // SublR is a rate [m/s]
  meltsub.st_energy.SublL=meltsub.st_energy.SublR*meltsub.st_energy.Dt;                                  // SublL=SublR*Dt; // SublL is a swe lenght [m]
}


//  FUNCTION TO CALCULATE MELT & SUBLIMATION RATE WHEN *NO* SNOW IS PRESSENT
void SurfaceEnergyBalance::EvapCalc (LocalData& evap) {
  evap.st_energy.EvR=evap.st_energy.fQe/(evap.st_energy.density_w*evap.st_energy.Le); // EvR=Qe/(ROWw*Le); // SublR is a rate [m/s]
  evap.st_energy.EvL=evap.st_energy.EvR*evap.st_energy.Dt;      // EvL=EvR*Dt;           // SublL is a swe lenght [m]
  if (evap.st_energy.ht_snow<=0.0000000000) {// Makesure Ml is not updated when a tiny amount of snow is present, Already taken care of by MeltSublRateCalc
    evap.st_energy.Ml = evap.st_energy.Pr+evap.st_energy.EvL;  // Ml = Pr+EvL;
    evap.st_energy.Mr=evap.st_energy.Ml/evap.st_energy.Dt;     // Mr=Ml/Dt;
  }
}


/*  ####     FUNCTION TO  CALCULATE THE CHANGE IN SNOWPACK DEPTH
    All Input variables are in (SWE) Snow Water equivialnce ***Even Ps (Precipitation as Snow) In the input files ***
    Here Ps gets converted to Snowpack dpeth and stays that way
    Zs, DELTz, TotwLoss is the Snowpack depth, ***NOT IN SWE EVER***
    To Convert from SWE to Snow pack depth (Zs) Multiply by density of water (ROWw) and Devide by Density of Snow (ROWs) or (ROWFs)
    SublL, Ml, Are always in SWE
*/
void SurfaceEnergyBalance::DeltaSnowPack (LocalData& deltsnow) {
  //if (SublL<0) {    //briefly changing SublL into snow (frost) lenght units *** This is for condensation ****
  //    CiL=-SublL*ROWw/ROWfrost, SublL=0 };  // I assume that ice is condensing on the snow surface therefore density is near that of ice
  //if (SublL>=0) {   //briefly changing SublL into snow lengh units *** This is for sublimation ****
  //    SublL=SublL*ROWw/ROWs };
  //TotwLoss = Ml*(ROWw/ROWs) + SublL - CiL;  //Calculate change in snow pack
  //if (SublL>=0) {  SublL=SublL*ROWs/ROWw };                     //Changing SublL back to SWE *** This is for sublimation ****
  if (deltsnow.st_energy.SublL<0) {    //briefly changing SublL into snow (frost) lenght units *** This is for condensation ****
    deltsnow.st_energy.CiL=-deltsnow.st_energy.SublL*deltsnow.st_energy.density_w/deltsnow.st_energy.density_frost;  // Assume ice is condensing on snow surface with a density near ice ~> 800[kg/m^3]
    deltsnow.st_energy.SublL=0;
  }
  if (deltsnow.st_energy.SublL>=0) {   //briefly changing SublL into snow lengh units *** This is for sublimation ****
    deltsnow.st_energy.SublL=deltsnow.st_energy.SublL*deltsnow.st_energy.density_w/deltsnow.st_energy.density_snow;
  }
  deltsnow.st_energy.TotwLoss = deltsnow.st_energy.Ml*(deltsnow.st_energy.density_w/deltsnow.st_energy.density_snow) + deltsnow.st_energy.SublL - deltsnow.st_energy.CiL;  //Calculate change in snow pack
  if (deltsnow.st_energy.SublL>=0) {                       //Changing SublL back to SWE *** This is for sublimation ****
    deltsnow.st_energy.SublL=deltsnow.st_energy.SublL*deltsnow.st_energy.density_freshsnow/deltsnow.st_energy.density_w;
  }
}


// FUNCTION TO MAKE SURE PROPER MASS OF WATER IS DELEVERED TO ATS WHEN SNOWPACK DISAPEARS
void SurfaceEnergyBalance::WaterMassCorr(LocalData& wmass) {
  //if (Ml>0) {
  //   if (SublL>0) {                             //When SublL is great then 0 Sublimation! Find the ratio of water atributed to SublL and subtract it from the water delivered to ATS
  //      SublL=(SublL*(ROWw/ROWs)/abs(TotwLoss));//SublL is temporarily changed to snow length
  //      SublL=SublL*Zs*(ROWs/ROWw);             //SublL is changed back to SWE & ratio of avaialbe snowpack for sublimation is found
  //      Ml=(Zs+Ps+(SublL*(ROWw/ROWs)))*(ROWs/ROWw);
  //   }else{Ml=(Zs+Ps+(SublL*-1))*(ROWs/ROWw) }; //SublL is less then 0, Frost! Add it to the water delivered to ATS
  //if (Ml==0) { SublL=(Zs+Ps)*(ROWs/ROWw) };     // The snowpack sublimated away.
  //TotwLoss=Zs;
  //Mr=Ml/Dt;}
  double Tswe = wmass.st_energy.density_w/wmass.st_energy.density_snow;
  double TZs = wmass.st_energy.density_snow/wmass.st_energy.density_w;
  if (wmass.st_energy.Ml>0) {
    if (wmass.st_energy.SublL>0) {               // When SublL is great then 0 Sublimation! Find the ratio of water atributed to SublL and subtract it from the water delivered to ATS

      wmass.st_energy.SublL=(wmass.st_energy.SublL*(Tswe)/std::abs(wmass.st_energy.TotwLoss));  //SublL is temporarily changed to snow length
      wmass.st_energy.SublL=wmass.st_energy.SublL*wmass.st_energy.ht_snow*(TZs);  //SublL is changed back to SWE & ratio of avaialbe snowpack for sublimation is found
      wmass.st_energy.Ml=(wmass.st_energy.ht_snow+wmass.st_energy.Ps+(wmass.st_energy.SublL*(Tswe)))*(TZs);
    } else {
      wmass.st_energy.Ml=(wmass.st_energy.ht_snow+wmass.st_energy.Ps+(wmass.st_energy.SublL*-1))*(TZs); //SublL is less then 0, Frost! Add it to water delivered to ATS
    }
  }
  if (wmass.st_energy.Ml==0) {
    wmass.st_energy.SublL=(wmass.st_energy.ht_snow+wmass.st_energy.Ps)*(TZs); // The snowpack sublimated away.
  }
  wmass.st_energy.TotwLoss=wmass.st_energy.ht_snow;
  wmass.st_energy.Mr=wmass.st_energy.Ml/wmass.st_energy.Dt;
}



/*   FUNCTION TO CALCULATE ENERGY BALANCE AND Ml WHEN SNOWPACK IS TEENY TINY
     This section is for to trim the enery deleiverd/taken from the soil with the snowpack is really really small
     Qc = -Ks*(Ts-Tb)/Zs;  --> blows up.  Boom!
     When Zs < 0.009 [m] or 1 cm we callcualte energy balance fromt he soil surface
     *****  This is ONLY done with the snow is melting ******
     **  Otherwise Snow Energy Balance is calucalted **
     */
void SurfaceEnergyBalance::TeenyTinySnowPack (LocalData& tiny) {
  tiny.st_energy.Qm=tiny.st_energy.Mr*(tiny.st_energy.density_w*tiny.st_energy.Hf); //Recalculating Qm based off of Snowpack limit
  double ZsHold = tiny.st_energy.ht_snow;
  tiny.st_energy.ht_snow=0.0;
  AlbedoCalc(tiny); // Calculating Surface Albedo (No Snow)
  tiny.st_energy.fQswIn=(1-tiny.st_energy.albedo_value)*tiny.st_energy.QswIn;
  GroundEnergyCalc(tiny); //  Solving Energy balance for the ground
  tiny.st_energy.funcall="BARE-teenytiny";
  tiny.st_energy.ht_snow=ZsHold;
  EvapCalc(tiny);//Calculating Evaporation from bare-ground "ht_snow" is keyed to *NOT* recalculate Ml because WaterMassCorr has already done that
  tiny.st_energy.Ml=tiny.st_energy.Ml+tiny.st_energy.Pr+tiny.st_energy.EvL;
}


// FUNCTION TO ADDS UP ALL THE CHANGES TO THE SNOWPACK
void SurfaceEnergyBalance::SnowPackCalc (LocalData& snowheight) {
  // Zs = Zs + Ps - TotwLoss;    // New formate
  snowheight.st_energy.ht_snow = snowheight.st_energy.ht_snow + snowheight.st_energy.Ps - snowheight.st_energy.TotwLoss;  //DELTz is the new snowpack depth
}


// FUNCTION TO TRACKS THE TIME (IN DAYS) WHERE NO NEW SNOW AS FALLEN ~> USED IN SNOW DENSITY
void SurfaceEnergyBalance::TrackSnowDays (LocalData& snowday) {
  if (snowday.st_energy.Ps<0.0001) {//If less then a mm of snow
    snowday.st_energy.nosnowdays=snowday.st_energy.nosnowdays+(snowday.st_energy.Dt/86400);
  }
  if (snowday.st_energy.Ps>=0.0001) {
    snowday.st_energy.nosnowdays=0;
  }
}



// FUNCTION TO CALCULATE SNOWPACK DENSITY ~> WEIGHTED AVERAGE OVER THREE POTEINTAL LAYERS OF SNOW
void SurfaceEnergyBalance::SnowPackDensity (LocalData& row) {
  if (row.st_energy.ht_snow<=0.0000000) {
    row.st_energy.density_snow=100;
    row.st_energy.HoldROWs=100;
  } else {
    // Track days with now snow and formulate snow deformation
    double ndensity = std::pow((row.st_energy.nosnowdays+1),0.3);
    row.st_energy.density_snow=row.st_energy.density_freshsnow*ndensity;
    //Weighted average of the three layers of the snowpack
    //ROWs=((Ps*ROWfs)+(Zs*ROWs)+(CiL*ROWfrost))/(Ps+Zs+CiL);  // Wt. average for snow density
    double denominator=row.st_energy.Ps+row.st_energy.ht_snow+row.st_energy.CiL;
    row.st_energy.density_snow=((row.st_energy.Ps*row.st_energy.density_freshsnow)+(row.st_energy.ht_snow*row.st_energy.density_snow)+(row.st_energy.CiL*row.st_energy.density_frost))/(denominator);
    if (row.st_energy.density_snow>950) {// Capping snow density ~> it should NEVER get this dense anyway.
      row.st_energy.density_snow=950.1;
    }
  }
}


// FUNCTION TO CONVERT TO SWE TO FRESHLY FALLEN SNOW DEPTH
void SurfaceEnergyBalance::SWE (LocalData& swe) {
  swe.st_energy.Ps = swe.st_energy.Ps*(swe.st_energy.density_w/swe.st_energy.density_freshsnow);
}

// FUNCTION TO FIND TEMPURATURE OF WATER FLOWING INTP ATS
void SurfaceEnergyBalance::WaterTemp(LocalData& wtemp) {
  if (wtemp.st_energy.Ml<=0.0) {
    wtemp.st_energy.Trw=wtemp.st_energy.air_temp;
  } else {
    wtemp.st_energy.Trw=273.15;
  }
}


// MAIN SNOW ENERGY BALANCE FUNCTION
void SurfaceEnergyBalance::SnowEnergyBalace (LocalData& seb) {
  // Calculated data/parameters
  double Qm=0.0, Ml=0.0;
  //Constant Parameters for Energy Equaions




  // Caculate Vapor pressure and dewpoint temperature from Air
  VaporCalc(seb.vp_air);
  //FILL IN DATA STRUCTURE NEEDED FOR PARAMETERS IN ENERGY BALANCE


  // Find Albedo
  AlbedoCalc(seb); // Calculating Surface Albedo
  seb.st_energy.fQswIn=(1-seb.st_energy.albedo_value)*seb.st_energy.QswIn;

  //Find Thermal Conductivity of snow
  ThermalConductSnow (seb); //Calculating themeral conductifity of snow

  //Convert snow leath Ps (m) in SWE lenght
  SWE(seb);

  if (seb.st_energy.ht_snow>0) {// ### IF THERE IS SNOW CHECK ### IF THERE IS SNOW CHECK ### IF THERE IS SNOW CHECK ##################
    BisectionEnergyCalc(seb); //Solving for snowsurface temperature and Energy balance equation
    seb.st_energy.funcall="SNOW";
    if (seb.st_energy.Ts<=273.15) {
      seb.st_energy.Qm=0;// Snow is not melting no water leaving snowpack as melt water
      seb.st_energy.Ml=0;
    }
    if (seb.st_energy.Ts>273.15) {// If snowsurface is warmer then ice
      seb.st_energy.Ts=273.15; // Setting snow temperature to zero
      MeltEnergyCalc(seb); // Recaculating Energy Balance when melt in occuring
      seb.st_energy.funcall="MELT";
    }//End if snowsurface is warmer then ice

    MeltSublRateCalc(seb); //Calculating melt, sublimation, and condisitation rate
    DeltaSnowPack(seb);  //CALCULATE THE CHANGE IN SNOWPACK DEPTH

    // Make sure proper mass of snowpack water gets delivered to ATS
    if (seb.st_energy.ht_snow < ((seb.st_energy.Ps - seb.st_energy.TotwLoss)*-1)) {
      WaterMassCorr(seb);
      seb.st_energy.funcall="Melt too small";
    }

    // CHECK FOR SMALL MELTING SNOWPACK, WHICH CAN CAUSE SUPPER HIGH QC Qc = -Ks*(Ts-Tb)/Zs;  --> blows up
    if (seb.st_energy.ht_snow<=0.009) {
      if (seb.st_energy.Ts==273.15) {
        TeenyTinySnowPack(seb);
      }
    }//  End section to trim energy deleivered when Qc = -Ks*(Ts-Tb)/Zs;  --> blows up

  }// ### END IF THERE IS SNOW CHECK  ### END IF THERE IS SNOW CHECK   ### END IF THERE IS SNOW CHECK ############################################

  if (seb.st_energy.ht_snow<=0) {// ### IF THERE IS NO SNOW CHECK  ###  ### IF THERE IS NO SNOW CHECK  ###  #################
    GroundEnergyCalc(seb);
    seb.st_energy.funcall="BARE";
    EvapCalc(seb);
  } // ### END IF THERE IS NO SNOW CHECK ###  ### END IF THERE IS NO SNOW CHECK ###  #######################

  //  UPDATE SNOWPACK DENSITY
  // NoSnowDays Updata
  TrackSnowDays(seb);
  // CALL SNOWPACK DENSITY FUNCTION
  seb.st_energy.HoldROWs=seb.st_energy.density_snow;
  SnowPackDensity(seb);

  // UDATE SNOWPACK HEIGHT
  seb.st_energy.HoldZs=seb.st_energy.ht_snow;
  SnowPackCalc(seb);
  if (seb.st_energy.ht_snow<=0.0000001) {
    seb.st_energy.ht_snow=0.0;
  }

  // CALCULATE WATER TEMPURATURE
  WaterTemp(seb);

} // END SEB Function ##############################################################################################################################
