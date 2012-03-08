#include "vanGenuchtenModel.hpp"
#include "math.h"
#include <iostream>

// vanGenuchtenModel::vanGenuchtenModel(int meshblock_, double m_, double alpha_, 
// 				     double sr_, double p_atm_) :
//   m(m_), alpha(alpha_), sr(sr_), p_atm(p_atm_)
// {
//   n = 1.0/(1.0-m);
//   set_mesh_block(meshblock_);
// }

vanGenuchtenModel::vanGenuchtenModel(std::string region_, double m_, double alpha_, 
				     double sr_, double p_atm_) :
  m(m_), alpha(alpha_), sr(sr_), p_atm(p_atm_)
{
  n = 1.0/(1.0-m);
  set_region(region_);
}
 

double vanGenuchtenModel::k_relative(double p)
{
  const double pc_transition ( 50.0 ); 

  double pc = p_atm - p; // capillary pressure

  if (pc > pc_transition) {
    double se = pow(1.0 + pow(alpha*pc,n),-m);
    
    return sqrt(se) * pow( 1.0 - pow( 1.0 - pow(se,1.0/m),m), 2);
  } else if (pc <= 0.0) {
    return 1.0;
  } else {
    
    double se_pc_transition(pow(1.0 + pow(alpha*pc_transition,n),-m));
    double f_pc_transition(sqrt(se_pc_transition) * pow( 1.0 - pow( 1.0 - pow(se_pc_transition,1.0/m),m), 2));
    
    return quadr_fn(pc, pc_transition, f_pc_transition);
  }

}


double vanGenuchtenModel::saturation(double p)
{
  const double pc_transition ( 50.0 ); 

  double pc = p_atm - p; // capillary pressure
  if (pc > pc_transition) {
    return pow(1.0 + pow(alpha*pc,n),-m) * (1.0 - sr) + sr;
  } else if (pc <= 0.0 ) {
    return 1.0;
  } else {
    
    double f_pc_transition(pow(1.0 + pow(alpha*pc_transition,n),-m) * (1.0 - sr) + sr);
    return quadr_fn(pc, pc_transition, f_pc_transition);
  }
}


double vanGenuchtenModel::d_saturation(double p)
{

  double pc = p_atm - p; // capillary pressure
  if (pc > 0.0) {
    return m * n * pow(1.0 + pow(alpha*pc,n),-m-1.0) * pow(alpha*pc,n-1.0) * alpha  * (1.0 - sr);
  } else { 
    return 0.0;
  } 
}  


double vanGenuchtenModel::pressure(double sl)
{
  double se = (sl - sr)/(1.0 - sr);
  return p_atm - ( pow( pow(se,-1.0/m) - 1.0, 1.0/n) )/alpha;
}


void vanGenuchtenModel::update_p_atm(double new_p_atm)
{
  p_atm = new_p_atm;
}


double vanGenuchtenModel::quadr_fn(double x, double pc_t, double f_pc_t) {
  
  return (f_pc_t - 1.0)/(pc_t*pc_t) * x*x + 1.0;

}
