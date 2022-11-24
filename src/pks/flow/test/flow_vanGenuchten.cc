/*
  Flow PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/

#include <iostream>
#include <cstdio>
#include <cmath>

#include "UnitTest++.h"

#include "WRM_vanGenuchten.hh"


TEST(vanGenuchten)
{
  using namespace Amanzi::Flow;

  double m = 0.5;
  double l = 0.5;
  double alpha = 0.01;
  double sr = 0.4;
  double p_atm = 101325.0;
  std::string krel_function("Mualem");
  double pc0 = 500.0;

  WRM_vanGenuchten vG(m, l, alpha, sr, krel_function, pc0);

  // check k_relative for p = 2*p_atm
  double pc = -p_atm;
  CHECK_EQUAL(vG.k_relative(pc), 1.0);

  // check k_relative for p = 0
  pc = p_atm;
  double se = pow(1.0 + pow(alpha * pc, 1.0 / (1.0 - m)), -m);
  CHECK_CLOSE(vG.k_relative(pc), sqrt(se) * pow(1.0 - pow(1.0 - pow(se, 1.0 / m), m), 2.0), 1e-15);

  // check saturation for p = 2*p_atm
  pc = -p_atm;
  CHECK_EQUAL(vG.saturation(pc), 1.0);

  // check saturation for p = 0
  pc = p_atm;
  CHECK_CLOSE(
    vG.saturation(pc), pow(1.0 + pow(alpha * pc, 1.0 / (1.0 - m)), -m) * (1.0 - sr) + sr, 1e-15);

  // check derivative of saturation(pc) at p = 2*p_atm.
  pc = -p_atm;
  CHECK_EQUAL(vG.dSdPc(pc), 0.0);

  // check derivative of saturation(p) at p = 0.
  pc = p_atm;
  CHECK_CLOSE(vG.dSdPc(pc),
              (1.0 - sr) * m * pow(1.0 + pow(alpha * pc, 1.0 / (1.0 - m)), -m - 1.0) * (-alpha) *
                pow(alpha * pc, m / (1.0 - m)) / (1.0 - m),
              1e-15);

  // check smoothing at p = 0.998 * p_atm
  pc = 0.002 * p_atm;
  CHECK_CLOSE(vG.k_relative(pc), 6.404827195147e-1, 1e-9);
}
