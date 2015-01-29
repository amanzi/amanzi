/*
Author: Ethan Coon

Painter's permafrost model with freezing point depression.

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
  double T0_ = 273.15; 
  double delT_ = 1.0;  
  double gamma_ = 7.33589e5; 
}

// required methods from the base class
// sats[0] = sg, sats[1] = sl, sats[2] = si
void
WRMFPDSmoothedPermafrostModel::saturations(double pc_liq, double pc_ice,
        double T, double (&sats)[3]) {

  double Tf = T0_*(1.-pc_liq/gamma_); //freezing temperature  
  double sl_sm =  T< Tf ? std::exp((T-Tf)/delT_) : 1. ;  
  sl_sm = ( wrm_->saturation(pc_liq) - wrm_->residualSaturation() )*sl_sm + wrm_->residualSaturation() ; 

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
  double Tf = T0_*(1.-pc_liq/gamma_); //freezing temperature  
  double f_ofT =  T< Tf ? std::exp((T-Tf)/delT_) : 1. ;  
  double sl_sm = ( wrm_->saturation(pc_liq) - wrm_->residualSaturation() )*f_ofT + wrm_->residualSaturation() ; 
  double sl = wrm_->saturation( std::max(pc_liq,pc_ice) ); 
  double sl_uf = wrm_->saturation( pc_liq ); //unfrozen conditions 
  if (pc_liq <= 0.) { // saturated
    dsats[0] = 0.; // gas
    if ( sl < sl_sm ) { 
      dsats[1] = wrm_->d_saturation(pc_liq)*f_ofT; 
      dsats[2] = -dsats[1];  
    } else { 
      dsats[1] = 0.;
      dsats[2] = 0.;
    } 
  } else if (pc_ice <= pc_liq) {
    dsats[2] = 0.; // ice
    dsats[1] = wrm_->d_saturation(pc_liq);  // liquid
    dsats[0] = - dsats[1];  // gas
  } else {
    if( sl < sl_sm ) {  
      dsats[1] = wrm_->d_saturation(pc_liq)*f_ofT; 
      dsats[2] = - dsats[1] *(1./sl_uf - sl_sm/std::pow(sl_uf,2));  
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
  double Tf = T0_*(1.-pc_liq/gamma_); //freezing temperature  
  double f_ofT =  T< Tf ? std::exp((T-Tf)/delT_) : 1. ;  
  double sl_sm = ( wrm_->saturation(pc_liq) - wrm_->residualSaturation() )*f_ofT + wrm_->residualSaturation() ; 
  double sl = wrm_->saturation( std::max(pc_liq,pc_ice) ); 
  double sl_uf = wrm_->saturation( pc_liq ); //unfrozen conditions 
  if (pc_liq <= 0.) { // saturated
    dsats[0] = 0.; // gas
    if ( sl < sl_sm ) { 
      dsats[1] = 0.; 
      dsats[2] = 0.;  
    } else { 
      dsats[1] = wrm_->d_saturation(pc_ice); // liquid
      dsats[2] = - dsats[1];  // ice
    } 
  } else if (pc_ice <= pc_liq) {
    dsats[2] = 0.; // ice
    dsats[1] = 0.;
    dsats[0] = 0.;
  } else {
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
  double Tf = T0_*(1.-pc_liq/gamma_); //freezing temperature  
  double f_ofT =  T< Tf ? std::exp((T-Tf)/delT_)/delT_ : 0. ;  
  double sl_sm = ( wrm_->saturation(pc_liq) - wrm_->residualSaturation() )*f_ofT + wrm_->residualSaturation() ; 
  double sl = wrm_->saturation( std::max(pc_liq,pc_ice) ); 
  double sl_uf = wrm_->saturation( pc_liq ); //unfrozen conditions 
  double sr = wrm_->residualSaturation(); 
  if ( sl > sl_sm || pc_liq > pc_ice ) { 
    dsats[0] = 0.;
    dsats[1] = 0.;
    dsats[2] = 0.;
  } else { 
    if( pc_liq <= 0. ) { 
      dsats[1] = (sl_uf-sr)*f_ofT;  
      dsats[2] = -dsats[1]; 
      dsats[0] = 0.; 
    } else { 
      dsats[1] = (sl_uf-sr)*f_ofT;  
      dsats[2] = -f_ofT/sl_uf; 
      dsats[0] = -dsats[1]-dsats[2]; 
    } 
  } 
 
}


} // namespace FlowRelations
} // namespace Flow
} // namespace Amanzi
