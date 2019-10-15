/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Markus Berndt
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

// TPLs
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "IOEvent.hh"


SUITE(IOEVENT)
{
  TEST(IOEVENT_CYCLES)
  {
    Teuchos::ParameterList plist;

    Teuchos::Array<int> p_cycles(std::vector<int>{ 0, 3, 4, 11 });
    plist.set("cycles", p_cycles);

    Amanzi::IOEvent V(plist);

    // test the cycle stuff, the expected result is in cycles_ and
    // we store the computed result in cycles
    int cycles_[13] = { 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0 };
    int cycles[13];
    for (int ic = 0; ic < 13; ic++) { cycles[ic] = V.DumpRequested(ic); }
    CHECK_ARRAY_EQUAL(cycles_, cycles, 13);
  }

  TEST(IOEVENT_CYCLES_SPS)
  {
    Teuchos::ParameterList plist;

    Teuchos::Array<int> csps(3);
    csps[0] = 0;
    csps[1] = 4;
    csps[2] = 10;
    plist.set<Teuchos::Array<int>>("cycles start period stop", csps);

    Amanzi::IOEvent V(plist);

    // test the cycle stuff, the expected result is in cycles_ and
    // we store the computed result in cycles
    int cycles_[31] = { 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    int cycles[31];
    for (int ic = 0; ic <= 30; ic++) { cycles[ic] = V.DumpRequested(ic); }
    CHECK_ARRAY_EQUAL(cycles_, cycles, 31);
  }


  TEST(IOEVENT_TIMES)
  {
    Teuchos::ParameterList plist;

    Teuchos::Array<double> tsps(3);
    tsps[0] = 0.0;
    tsps[1] = 4.0;
    tsps[2] = 10.0;
    plist.set<Teuchos::Array<double>>("times start period stop", tsps);

    Teuchos::Array<double> times(2);
    times[0] = 1.0;
    times[1] = 3.0;
    plist.set<Teuchos::Array<double>>("times", times);

    Amanzi::IOEvent V(plist);

    // test the time sps stuff
    CHECK_EQUAL(true, V.DumpRequested(0.0));
    CHECK_EQUAL(true, V.DumpRequested(1.0));
    CHECK_EQUAL(true, V.DumpRequested(3.0));
    CHECK_EQUAL(true, V.DumpRequested(4.0));
    CHECK_EQUAL(true, V.DumpRequested(8.0));

    CHECK_EQUAL(false, V.DumpRequested(0.5));
    CHECK_EQUAL(false, V.DumpRequested(1.1));
    CHECK_EQUAL(false, V.DumpRequested(3.2));
    CHECK_EQUAL(false, V.DumpRequested(3.99));
    CHECK_EQUAL(false, V.DumpRequested(10.0));
    CHECK_EQUAL(false, V.DumpRequested(12.0));
  }

  TEST(IOEVENT_MULTIPLE)
  {
    Teuchos::ParameterList plist;

    Teuchos::Array<double> tsps(3);
    tsps[0] = 0.0;
    tsps[1] = 4.0;
    tsps[2] = 10.0;
    plist.set<Teuchos::Array<double>>("times start period stop", tsps);

    Teuchos::Array<int> p_cycles(std::vector<int>{ 0, 3, 4, 11 });
    plist.set("cycles", p_cycles);

    Amanzi::IOEvent V(plist);

    // test the time sps stuff
    CHECK_EQUAL(true, V.DumpRequested(0.0));
    CHECK_EQUAL(true, V.DumpRequested(4.0));
    CHECK_EQUAL(true, V.DumpRequested(8.0));

    CHECK_EQUAL(false, V.DumpRequested(0.5));
    CHECK_EQUAL(false, V.DumpRequested(1.1));
    CHECK_EQUAL(false, V.DumpRequested(3.2));
    CHECK_EQUAL(false, V.DumpRequested(3.99));
    CHECK_EQUAL(false, V.DumpRequested(10.0));
    CHECK_EQUAL(false, V.DumpRequested(12.0));


    // test the cycle stuff, the expected result is in cycles_ and
    // we store the computed result in cycles
    int cycles_[13] = { 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0 };
    int cycles[13];
    for (int ic = 0; ic < 13; ic++) { cycles[ic] = V.DumpRequested(ic); }
    CHECK_ARRAY_EQUAL(cycles_, cycles, 13);

    // test the two-parameter call
    CHECK_EQUAL(true, V.DumpRequested(1, 0.0));
    CHECK_EQUAL(true, V.DumpRequested(1, 4.0));
    CHECK_EQUAL(true, V.DumpRequested(1, 8.0));
    CHECK_EQUAL(false, V.DumpRequested(1, 0.5));
    CHECK_EQUAL(false, V.DumpRequested(1, 1.1));
    CHECK_EQUAL(false, V.DumpRequested(1, 3.2));
    CHECK_EQUAL(false, V.DumpRequested(1, 3.99));
    CHECK_EQUAL(false, V.DumpRequested(1, 10.0));
    CHECK_EQUAL(false, V.DumpRequested(1, 12.0));

    for (int ic = 0; ic < 13; ic++) { cycles[ic] = V.DumpRequested(ic, 11.0); }
    CHECK_ARRAY_EQUAL(cycles_, cycles, 13);

    CHECK_EQUAL(true, V.DumpRequested(0, 0.0));
    CHECK_EQUAL(true, V.DumpRequested(3, 4.0));
  }

  TEST(IOEVENT_TIMES_WITH_UNITS)
  {
    Teuchos::ParameterList plist;

    Teuchos::Array<double> tsps(3);
    tsps[0] = 0.0;
    tsps[1] = 4.0;
    tsps[2] = 10.0;
    plist.set<Teuchos::Array<double>>("times start period stop", tsps);
    plist.set<std::string>("times start period stop units", "noleap");

    Teuchos::Array<double> times(2);
    times[0] = 1.0;
    times[1] = 3.0;
    plist.set<Teuchos::Array<double>>("times", times);
    plist.set<std::string>("times units", "d");

    Amanzi::IOEvent V(plist);

    // test the time sps stuff
    double y_s = 365. * 86400.;
    double d_s = 86400.;

    CHECK_EQUAL(true, V.DumpRequested(0.0));
    CHECK_EQUAL(true, V.DumpRequested(1.0 * d_s));
    CHECK_EQUAL(true, V.DumpRequested(3.0 * d_s));
    CHECK_EQUAL(true, V.DumpRequested(4.0 * y_s));
    CHECK_EQUAL(true, V.DumpRequested(8.0 * y_s));

    CHECK_EQUAL(false, V.DumpRequested(0.5 * d_s));
    CHECK_EQUAL(false, V.DumpRequested(1.1 * d_s));
    CHECK_EQUAL(false, V.DumpRequested(3.2 * d_s));
    CHECK_EQUAL(false, V.DumpRequested(3.99 * y_s));
    CHECK_EQUAL(false, V.DumpRequested(10.0 * y_s));
    CHECK_EQUAL(false, V.DumpRequested(12.0 * y_s));
  }


  TEST(IOEVENT_TIMES_WITH_UNITS_NEGATIVE)
  {
    Teuchos::ParameterList plist;

    Teuchos::Array<double> tsps(3);
    tsps[0] = 0.0;
    tsps[1] = 4.0;
    tsps[2] = -1.0;
    plist.set<Teuchos::Array<double>>("times start period stop", tsps);
    plist.set<std::string>("times start period stop units", "noleap");

    Teuchos::Array<double> times(2);
    times[0] = 1.0;
    times[1] = 3.0;
    plist.set<Teuchos::Array<double>>("times", times);
    plist.set<std::string>("times units", "d");

    Amanzi::IOEvent V(plist);

    // test the time sps stuff
    double y_s = 365. * 86400.;
    double d_s = 86400.;

    CHECK_EQUAL(true, V.DumpRequested(0.0));
    CHECK_EQUAL(true, V.DumpRequested(1.0 * d_s));
    CHECK_EQUAL(true, V.DumpRequested(3.0 * d_s));
    CHECK_EQUAL(true, V.DumpRequested(4.0 * y_s));
    CHECK_EQUAL(true, V.DumpRequested(8.0 * y_s));
    CHECK_EQUAL(true, V.DumpRequested(12.0 * y_s));

    CHECK_EQUAL(false, V.DumpRequested(0.5 * d_s));
    CHECK_EQUAL(false, V.DumpRequested(1.1 * d_s));
    CHECK_EQUAL(false, V.DumpRequested(3.2 * d_s));
    CHECK_EQUAL(false, V.DumpRequested(3.99 * y_s));
    CHECK_EQUAL(false, V.DumpRequested(10.0 * y_s));
  }
}
