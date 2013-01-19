#include <iostream>
#include "UnitTest++.h"

#include "wrm_van_genuchten.hh"
#include "math.h"

TEST(vanGenuchten) {
  using namespace Amanzi::Flow::FlowRelations;

  double m = 0.5;
  double alpha = 0.1;
  double sr = 0.4;
  double p_atm = 1.0e+5;

  Teuchos::ParameterList plist;
  plist.set("van Genuchten m", m);
  plist.set("van Genuchten alpha", alpha);
  plist.set("van Genuchten residual saturation", sr);
  plist.set("van Genuchten smoothing interval width", 0.0);
  WRMVanGenuchten vG(plist);

  // check k_relative for p = 2*p_atm
  double pc = -p_atm;
  CHECK_EQUAL(vG.k_relative(pc), 1.0);

  // check k_relative for p = 0
  pc = p_atm;
  double se = pow(1.0 + pow(alpha * pc, 1.0 / (1.0-m)),-m);
  CHECK_CLOSE(vG.k_relative(pc),
              sqrt(se) * pow(1.0 - pow(1.0 - pow(se, 1.0/m), m), 2.0), 1e-15);

  // check saturation for p = 2*p_atm
  pc = -p_atm;
  CHECK_EQUAL(vG.saturation(pc), 1.0);

  // check saturation for p = 0
  pc = p_atm;
  CHECK_CLOSE(vG.saturation(pc),
              pow(1.0 + pow(alpha * pc, 1.0/ (1.0-m)), -m) * (1.0-sr) + sr, 1e-15);

  // check derivative of saturation(p) at p=2*p_atm
  pc = -p_atm;
  CHECK_EQUAL(vG.d_saturation(pc), 0.0);

  // check derivative of saturation(p) at p=0
  pc = p_atm;
  CHECK_CLOSE(vG.d_saturation(pc),
              (1.0-sr)*(-m)*pow(1.0+pow(alpha*pc,1.0/(1.0-m)),-m-1.0)
              *alpha*pow(alpha*pc,m/(1.0-m))/(1.0-m),1e-15);

  // check capillary pressure at p = 2*p_atm
  pc = -p_atm;
  CHECK_CLOSE(vG.capillaryPressure( vG.saturation(pc) ), pc, 1e-15);

  // check capillary pressure at p = 0.
  pc = p_atm;
  CHECK_CLOSE(vG.capillaryPressure( vG.saturation(pc) ), pc, 1e-15);

  // check d_capillaryPressure at p = 2*p_atm
  pc = -p_atm;
  CHECK_CLOSE(vG.d_capillaryPressure( vG.saturation(pc) ),
              1.0 / vG.d_saturation(pc), 1e-15);

  // check d_capillaryPressure at p = 0
  pc = p_atm;
  CHECK_CLOSE(vG.d_capillaryPressure( vG.saturation(pc) ),
              1.0 / vG.d_saturation(pc), 1e-15);
}
