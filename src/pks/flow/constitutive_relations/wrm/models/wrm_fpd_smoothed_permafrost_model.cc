/*
Author: Ethan Coon, Scott Painter 

Painter's permafrost model with freezing point depression.
Smoothed to make life easier. 

 */

#include "dbc.hh"

#include "wrm.hh"
#include "wrm_fpd_smoothed_permafrost_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// Constructor
WRMFPDSmoothedPermafrostModel::WRMFPDSmoothedPermafrostModel(Teuchos::ParameterList& plist) :
     WRMPermafrostModel(plist) {

  T0_ = plist_.get<double>("heat of fusion reference temperature [K]", 273.15);
  delT_ = plist_.get<double>("width of smoothing interval [K]", 1.0);

  double sigma_ice_liq = plist_.get<double>("interfacial tension ice-water", 33.1);
  double sigma_gas_liq = plist_.get<double>("interfacial tension air-water", 72.7);
  double heat_fusion = plist_.get<double>("heat of fusion of water [J/kg]", 3.34e5);
  double dens = 999.915; // water density [kg/m3] at 0 C

  gamma_ = sigma_gas_liq/sigma_ice_liq * heat_fusion * dens; 
}


// required methods from the base class
// sats[0] = sg, sats[1] = sl, sats[2] = si
void
WRMFPDSmoothedPermafrostModel::saturations(double pc_liq, double pc_ice,
        double T, double (&sats)[3]) {

  double Tf = T0_*(1.-std::max(0., pc_liq)/gamma_); //freezing temperature  
  double sr = wrm_->residualSaturation(); // residiual  
//partitioning function in transition zone 
  double sl_sm =  T< Tf ? std::exp((T-Tf)/delT_) : 1. ;  
  sl_sm = ( wrm_->saturation(pc_liq) - sr )*sl_sm + sr; 

  if (pc_liq <= 0.) { // saturated
    sats[1] = wrm_->saturation(pc_ice); // liquid
    if( sats[1] < sl_sm ) sats[1]=sl_sm; 
    sats[2] = 1.0 - sats[1];  // ice
    sats[0] = 0.; // gas
    ASSERT(sats[2] >= 0.);
  } else if (pc_ice <= pc_liq) {
    sats[1] = wrm_->saturation(pc_liq);  // liquid
    sats[2] = 0.; // ice
    sats[0] = 1.0 - sats[1];  // gas
    ASSERT(sats[0] >= 0.);
  } else {
    sats[1] = wrm_->saturation(pc_ice);
    if( sats[1] < sl_sm ) sats[1]=sl_sm;  
    sats[2] = 1. - sats[1] / wrm_->saturation(pc_liq);
    sats[0] = 1. - sats[1] - sats[2];
    ASSERT(sats[2] >= 0.);
    ASSERT(sats[0] >= 0.);
  }

    
}

void
WRMFPDSmoothedPermafrostModel::dsaturations_dpc_liq(double pc_liq, double pc_ice,
        double T, double (&dsats)[3]) {
  double sl = 0.; 
  double sl_uf = wrm_->saturation( pc_liq ); //unfrozen conditions 
  double Tf = T0_*(1.-std::max(0., pc_liq)/gamma_); //freezing temperature  
  double f_ofT =  T< Tf ? std::exp((T-Tf)/delT_) : 1. ;  
  double sr = wrm_->residualSaturation(); 
  double sl_sm = ( sl_uf - sr )*f_ofT + sr; 
  if (pc_liq <= 0.) { // saturated
    sl = wrm_->saturation(pc_ice); 
    dsats[0] = 0.; // gas
    if ( sl < sl_sm ) { //transition zone  
      dsats[1] = wrm_->d_saturation(pc_liq)*f_ofT + sl_uf*T0_*f_ofT/(gamma_*delT_); 
      dsats[2] = -dsats[1];  
    } else { 
      dsats[1] = 0.;
      dsats[2] = 0.;
    } 
  } else if (pc_ice <= pc_liq) {
    sl = wrm_->saturation(pc_liq); 
    dsats[2] = 0.; // ice
    dsats[1] = wrm_->d_saturation(pc_liq);  // liquid
    dsats[0] = - dsats[1];  // gas
  } else {
    sl = wrm_->saturation(pc_ice); 
    if( sl < sl_sm ) { //transition zone  
      dsats[1] = wrm_->d_saturation(pc_liq)*f_ofT + sl_uf*T0_*f_ofT/(gamma_*delT_); 
      dsats[2] = (sl_sm/sl_uf*wrm_->d_saturation(pc_liq) - dsats[1])/sl_uf; 
      dsats[0] = -dsats[1] - dsats[2]; 
    } else { 
      dsats[1] = 0.;
      dsats[2] = wrm_->saturation(pc_ice) / std::pow(wrm_->saturation(pc_liq),2) 
         * wrm_->d_saturation(pc_liq);
      dsats[0] = - dsats[2];
    }      
  }
}


void
WRMFPDSmoothedPermafrostModel::dsaturations_dpc_ice(double pc_liq, double pc_ice,
        double T, double (&dsats)[3]) {
  double sl = 0.;  
  double sl_uf = wrm_->saturation( pc_liq ); //unfrozen conditions 
  double Tf = T0_*(1.-std::max(0., pc_liq)/gamma_); //freezing temperature  
  double f_ofT =  T< Tf ? std::exp((T-Tf)/delT_) : 1. ;  
  double sr = wrm_->residualSaturation(); 
  double sl_sm = ( sl_uf - sr )*f_ofT + sr; 
  if (pc_liq <= 0.) { // saturated
    dsats[0] = 0.; // gas
    sl = wrm_->saturation(pc_ice); 
    if ( sl < sl_sm ) { 
      dsats[1] = 0.; 
      dsats[2] = 0.;  
    } else { 
      dsats[1] = wrm_->d_saturation(pc_ice); // liquid
      dsats[2] = - dsats[1];  // ice
    } 
  } else if (pc_ice <= pc_liq) {
    sl = wrm_->saturation(pc_liq); 
    dsats[2] = 0.; // ice
    dsats[1] = 0.;
    dsats[0] = 0.;
  } else {
    sl = wrm_->saturation(pc_ice); 
    if( sl < sl_sm) { 
      dsats[1] = 0.;  
      dsats[2] = 0.;  
      dsats[0] = 0.;  
    } else { 
      dsats[1] = wrm_->d_saturation(pc_ice);
      dsats[2] = - dsats[1] / wrm_->saturation(pc_liq);
      dsats[0] = - dsats[1] - dsats[2];
    } 
  }


}

void
WRMFPDSmoothedPermafrostModel::dsaturations_dtemperature(double pc_liq, double pc_ice,
        double T, double (&dsats)[3]) {
  double sl = wrm_->saturation( std::max(pc_liq,pc_ice) ); //frozen or unfrozen
  double sl_uf = wrm_->saturation( pc_liq ); //unfrozen conditions 
  double Tf = T0_*(1.-std::max(0., pc_liq)/gamma_); //freezing temperature  
  double f_ofT =  T< Tf ? std::exp((T-Tf)/delT_) : 1. ;  
  double sr = wrm_->residualSaturation(); 
  double sl_sm = ( sl_uf - sr )*f_ofT + sr; 
  double df_ofT =  T< Tf ? std::exp((T-Tf)/delT_)/delT_ : 0. ;  

  if (pc_liq <= 0.) { // saturated
    sl = wrm_->saturation(pc_ice); 
    if ( sl < sl_sm) { 
      dsats[1] = (sl_uf-sr)*df_ofT;  
      dsats[2] = -dsats[1]; 
      dsats[0] = 0.; 
    } else { 
      dsats[0] = 0.;
      dsats[1] = 0.;
      dsats[2] = 0.;
    } 
  } else if (pc_ice <= pc_liq) {
    dsats[0] = 0.;
    dsats[1] = 0.;
    dsats[2] = 0.;
  } else {
    sl = wrm_->saturation(pc_ice); 
    if ( sl < sl_sm) { 
      dsats[1] = (sl_uf-sr)*df_ofT;  
      dsats[2] = -dsats[1]/sl_uf; 
      dsats[0] = -dsats[1]-dsats[2]; 
    } 
  } 
 
}


} // namespace FlowRelations
} // namespace Flow
} // namespace Amanzi
