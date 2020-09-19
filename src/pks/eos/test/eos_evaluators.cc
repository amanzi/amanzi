/*
  EOS

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.
*/

#include <iostream>
#include <cstdio>
#include <cmath>

#include "UnitTest++.h"

#include "EOS_Water.hh"
#include "EOS_WaterFEHM.hh"


TEST(DensityEvaluators) {
  using namespace Amanzi::AmanziEOS;

  Teuchos::ParameterList plist;
  EOS_Water eos(plist);
 
  double p(101325.0), T0(273.15), d0, d1;
  for (double T = 0; T < 30; T += 0.1) {
    double d1 = eos.MassDensity(T + T0, p);
    CHECK_CLOSE(d1, 998.0, 2.0);

    if (T > 10.0) {
      CHECK(d1 < d0);
      double der = eos.DMassDensityDT(T + T0, p);
      CHECK(der < 0.0);
    }
    d0 = d1;
  }

  // next EOS
  EOS_WaterFEHM eos_fehm(plist);
 
  p = 101325.0;
  for (int loop = 0; loop < 5; loop++) {
    p *= 2.0;
    d0 = 1200.0;
    for (double T = 15; T < 300; T += 1.5) {
      double d1 = eos_fehm.MassDensity(T + T0, p);
      CHECK(d1 < 1001.0 && d1 > 550.0 && d1 < d0);
      d0 = d1;

      double der = eos_fehm.DMassDensityDT(T + T0, p);
      CHECK(der < 0.0);
    }
  }
}


