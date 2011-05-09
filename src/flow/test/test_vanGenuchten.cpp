#include <iostream>
#include "UnitTest++.h"

#include "vanGenuchtenModel.hpp"
#include "math.h"

TEST(vanGenuchten) {

  double m = 0.5;
  double alpha = 0.1;
  double sr = 0.4;
  double p_atm = 1.0e+5;

  vanGenuchtenModel vG(m, alpha, sr, p_atm);
  
  // check k_relative for p = 2*p_atm
  CHECK_EQUAL(vG.k_relative(2.0*p_atm),1.0);
  
  // check k_relative for p = 0, then pc = -p_atm
  double se = pow(1.0 + pow(alpha*p_atm,1.0/(1.0-m)),-m);
  CHECK_CLOSE(vG.k_relative(0.0),
	      sqrt(se)*pow(1.0-pow(1.0-pow(se,1.0/m),m),2.0),1e-15);
  
  // check saturation for p = 2*p_atm
  CHECK_EQUAL(vG.saturation(2*p_atm),1.0);
  
  // check saturation for p = 0, then pc = -p_atm
  CHECK_CLOSE(vG.saturation(0.0),
	      pow(1.0 + pow(alpha*p_atm,1.0/(1.0-m)),-m)*(1.0-sr)+sr,1e-15);
 

  // check derivative of saturation(p) at p=2*p_atm
  CHECK_EQUAL(vG.d_saturation(2*p_atm), 0.0);

  // check derivative of saturation(p) at p=0.0
  CHECK_CLOSE(vG.d_saturation(0.0), 
	      (1.0-sr)*(-m)*pow(1.0+pow(alpha*p_atm,1.0/(1.0-m)),-m-1.0)
	      *(-alpha)*pow(alpha*p_atm,m/(1.0-m))/(1.0-m),1e-15);
  
}


