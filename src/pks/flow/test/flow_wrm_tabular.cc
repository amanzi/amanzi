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

#include "WRM_tabular.hh"


TEST(WRM_TABULAR)
{
  using namespace Amanzi::Flow;

  Teuchos::ParameterList plist;
  plist.set<double>("tolerance", 1.0e-14);
  plist.set<Teuchos::Array<double>>("cap pressure",
                                    std::vector<double>({ 0.0, 100.0, 500.0, 1000.0, 10000.0 }));
  plist.set<Teuchos::Array<double>>("permeability",
                                    std::vector<double>({ 1.0, 0.5, 0.1, 0.05, 0.02 }));
  plist.set<Teuchos::Array<double>>("saturation",
                                    std::vector<double>({ 1.0, 0.5, 0.1, 0.05, 0.02 }));

  WRM_tabular wrm(plist);

  CHECK_EQUAL(wrm.k_relative(0.0), 1.0);

  for (double pc = 0.0; pc < 10000.0; pc += 123.0) {
    double s = wrm.saturation(pc);
    double pc_new = wrm.capillaryPressure(s);
    CHECK_CLOSE(pc, pc_new, 1e-10 * (pc + 1));
    CHECK(wrm.get_itrs() < 30);
  }
}
