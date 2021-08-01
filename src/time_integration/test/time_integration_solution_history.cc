/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "AmanziTypes.hh"
#include "AmanziComm.hh"
#include "AmanziVector.hh"
#include "SolutionHistory.hh"

using namespace Amanzi;

SUITE(SolutionHistoryTests)
{
  // data structures for testing
  struct test_data {
    Comm_ptr_type comm;
    Vector_ptr_type x;

    test_data()
    {
      comm = getDefaultComm();
      auto map = Teuchos::rcp(new Map_type(2, 0, comm));
      x = Teuchos::rcp(new Vector_type(map));
    }
  };

  TEST_FIXTURE(test_data, SolutionHistory_1)
  {
    std::cout << "Test: SolutionHistory_1" << std::endl;

    // create a solution history of size three
    Teuchos::RCP<Amanzi::SolutionHistory<Vector_type>> SH =
      Teuchos::rcp(new Amanzi::SolutionHistory<Vector_type>(3, 0.0, *x));
    x->putScalar(1.0);
    CHECK_EQUAL(SH->history_size(), 1);
    CHECK_EQUAL(SH->MostRecentTime(), 0.0);

    SH->RecordSolution(1.0, *x);

    CHECK_EQUAL(SH->history_size(), 2);
    CHECK_EQUAL(SH->MostRecentTime(), 1.0);

    x->putScalar(4.0);
    SH->RecordSolution(2.0, *x);

    CHECK_EQUAL(SH->history_size(), 3);
    CHECK_EQUAL(SH->MostRecentTime(), 2.0);

    x->putScalar(9.0);
    SH->RecordSolution(3.0, *x);

    CHECK_EQUAL(SH->history_size(), 3);
    CHECK_EQUAL(SH->MostRecentTime(), 3.0);

    // check that the most recent vector in fact is the
    // one that contains  9.0
    Teuchos::RCP<Vector_type> y = Teuchos::rcp(new Vector_type(x->getMap()));
    SH->MostRecentSolution(*y);

    double norminf =
      y->normInf(); // compute the max norm of the computed difference
    CHECK_EQUAL(norminf, 9.0);

    std::vector<double> h;
    SH->TimeDeltas(h);

    CHECK_EQUAL(h.size(), SH->history_size() - 1);
    CHECK_EQUAL(h[0], 1.0);
    CHECK_EQUAL(h[1], 2.0);

    // interpolate 1st order
    SH->InterpolateSolution(2.5, *y, 1);

    // we should get 6.5
    norminf = y->normInf();
    CHECK_EQUAL(norminf, 6.5);

    // interpolate 2nd order
    SH->InterpolateSolution(2.5, *y, 2);

    // from the quadratic through (1,1),(2,4),(3,9)
    // we should get 6.25
    // x->putScalar(6.25);
    // x->update(-1.0, *y, 1.0);
    norminf = y->normInf();
    CHECK_EQUAL(norminf, 6.25);

    // interpolate maximum order (should be 2nd)
    SH->InterpolateSolution(2.5, *y);

    // from the quadratic through (1,1),(2,4),(3,9)
    // we should get 6.25
    norminf = y->normInf();
    CHECK_EQUAL(norminf, 6.25);
  }

  TEST_FIXTURE(test_data, SolutionHistory_2)
  {
    std::cout << "Test: SolutionHistory_2" << std::endl;
    auto xdot = Teuchos::rcp(new Vector_type(x->getMap()));
    x->putScalar(0.0);
    xdot->putScalar(0.0);
    // create a solution history of size three
    Teuchos::RCP<Amanzi::SolutionHistory<Vector_type>> SH =
      Teuchos::rcp(new Amanzi::SolutionHistory<Vector_type>(4, 0.0, *x, *xdot));

    x->putScalar(1.0);
    xdot->putScalar(2.0);
    SH->RecordSolution(1.0, *x, *xdot);

    auto y = Teuchos::rcp(new Vector_type(x->getMap()));

    // interpolate 1st order
    SH->InterpolateSolution(0.5, *y, 2);

    double norminf = y->normInf();
    CHECK_EQUAL(norminf, 0.25);

    // interpolate maximum order (3rd in this case)
    SH->InterpolateSolution(0.5, *y);

    norminf = y->normInf();
    CHECK_EQUAL(norminf, 0.25);
  }
}
