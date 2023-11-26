/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Flow PK

*/

#include <iostream>
#include <cstdio>
#include <cmath>

#include "UnitTest++.h"
#include "Teuchos_ParameterList.hpp"

#include "ApertureModel_BartonBandis.hh"
#include "ApertureModel_ExponentialLaw.hh"

TEST(vanGenuchten)
{
  using namespace Amanzi;
  using namespace Amanzi::Flow;

  Teuchos::ParameterList plist;
  plist.set<double>("undeformed aperture", 1e-5)
    .set<double>("overburden pressure", 11e+6)
    .set<double>("BartonBandis A", 1e-11)
    .set<double>("BartonBandis B", 0.0)
    .set<double>("compressibility", 2e-7);

  ApertureModel_BartonBandis abb(plist);
  ApertureModel_ExponentialLaw ael(plist);

  for (double p = 10e+6; p < 30e+6; p += 1e+6) {
    double a1 = abb.Aperture(p);
    double a2 = ael.Aperture(p);
    CHECK(a1 > 0.0 && a2 > 0.0);
    CHECK(a1 < 1e-3 && a2 < 1e-3);
    // std::cout << p << " " << a1 << " " << a2 << std::endl;
  }
}
