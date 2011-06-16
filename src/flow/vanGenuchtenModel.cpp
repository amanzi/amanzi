#include "vanGenuchtenModel.hpp"
#include "math.h"
#include <iostream>

vanGenuchtenModel::vanGenuchtenModel(int meshblock_, double m_, double alpha_, 
				     double sr_, double p_atm_) :
  m(m_), alpha(alpha_), sr(sr_), p_atm(p_atm_)
{
  n = 1.0/(1.0-m);
  set_mesh_block(meshblock_);
}

 
double vanGenuchtenModel::k_relative(double p)
{
  double pc = p - p_atm; // capillary pressure
  if (pc < 0.0)
    {
      double se = pow(1.0 + pow(-alpha*pc,n),-m);
      return sqrt(se) * pow( 1.0 - pow( 1.0 - pow(se,1.0/m),m), 2.0);
    }
  else
    {
      return 1.0;
    }
  
}

double vanGenuchtenModel::saturation(double p)
{
  double pc = p - p_atm; // capillary pressure
  if (pc < 0.0) 
    {
      return pow(1.0 + pow(-alpha*pc,n),-m) * (1.0 - sr) + sr;
    }
  else
    {
      return 1.0;
    }
}

double vanGenuchtenModel::d_saturation(double p)
{
  double pc = p - p_atm; // capillary pressure
  if (pc < 0.0) 
    {
      return -m * pow(1.0 + pow(-alpha*pc,n),-m-1.0) 
	* n * pow(-alpha*pc,n - 1.0) * (-alpha) * (1.0 - sr);
    }
  else
    {
      return 0.0;
    }
}  


double vanGenuchtenModel::pressure(double sl)
{
  double se = (sl - sr)/(1.0 - sr);
  
  return p_atm - ( pow( pow(se,-1.0/m) - 1.0, 1/n) )/alpha;
  
}


void vanGenuchtenModel::update_p_atm(double new_p_atm)
{
  p_atm = new_p_atm;
}
