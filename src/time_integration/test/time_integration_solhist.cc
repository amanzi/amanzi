/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Vector.h"

#include "SolutionHistory.hh"

using namespace Amanzi;
using SolutionHistory_type = Amanzi::SolutionHistory<Epetra_Vector, Epetra_Map>;

SUITE(SolutionHistoryTests)
{
  // data structures for testing
  struct test_data {
    Comm_ptr_type comm;
    Teuchos::RCP<Epetra_Map> map;
    Teuchos::RCP<Epetra_Vector> x;

    test_data()
    {
      comm = getDefaultComm();
      map = Teuchos::rcp(new Epetra_Map(2, 0, *comm));
      x = Teuchos::rcp(new Epetra_Vector(*map));
    }
  };

  TEST_FIXTURE(test_data, SolutionHistory_1)
  {
    std::cout << "Test: SolutionHistory_1" << std::endl;

    // create a solution history of size three
    auto SH = Teuchos::rcp(new SolutionHistory_type("myhist", 3, 0.0, map));
    x->PutScalar(1.0);
    SH->FlushHistory(0., *x);

    CHECK_EQUAL(SH->history_size(), 1);
    CHECK_EQUAL(SH->MostRecentTime(), 0.0);

    SH->RecordSolution(1.0, *x);

    CHECK_EQUAL(SH->history_size(), 2);
    CHECK_EQUAL(SH->MostRecentTime(), 1.0);

    x->PutScalar(4.0);
    SH->RecordSolution(2.0, *x);

    CHECK_EQUAL(SH->history_size(), 3);
    CHECK_EQUAL(SH->MostRecentTime(), 2.0);

    x->PutScalar(9.0);
    SH->RecordSolution(3.0, *x);

    CHECK_EQUAL(SH->history_size(), 3);
    CHECK_EQUAL(SH->MostRecentTime(), 3.0);

    // check that the most recent vector in fact is the
    // one that contains  9.0
    auto y = Teuchos::rcp(new Epetra_Vector(*x));
    *y = *SH->MostRecentSolution();
    double norminf[1];
    y->NormInf(norminf); // compute the max norm of the computed difference
    CHECK_EQUAL(9.0, *norminf);

    std::vector<double> h;
    SH->TimeDeltas(h);

    CHECK_EQUAL(h.size(), SH->history_size() - 1);
    CHECK_EQUAL(1.0, h[0]);
    CHECK_EQUAL(2.0, h[1]);

    // interpolate 1st order
    SH->InterpolateSolution(2.5, *y, 1);

    // we should get 6.5
    y->NormInf(norminf);
    CHECK_EQUAL(6.5, norminf[0]);

    // interpolate 2nd order
    SH->InterpolateSolution(2.5, *y, 2);

    // from the quadratic through (1,1),(2,4),(3,9)
    // we should get 6.25
    // x->PutScalar(6.25);
    // x->Update(-1.0, *y, 1.0);
    y->NormInf(norminf);
    CHECK_EQUAL(6.25, norminf[0]);

    // interpolate maximum order (should be 2nd)
    SH->InterpolateSolution(2.5, *y);

    // from the quadratic through (1,1),(2,4),(3,9)
    // we should get 6.25
    y->NormInf(norminf);
    CHECK_EQUAL(6.25, norminf[0]);
  }

  TEST_FIXTURE(test_data, SolutionHistory_2)
  {
    std::cout << "Test: SolutionHistory_2" << std::endl;
    Teuchos::RCP<Epetra_Vector> xdot = Teuchos::rcp(new Epetra_Vector(*x));

    x->PutScalar(0.0);
    xdot->PutScalar(0.0);
    // create a solution history of size three
    auto SH = Teuchos::rcp(new SolutionHistory_type("myhist", 4, 0.0, map));
    SH->FlushHistory(0., *x, xdot.get());

    x->PutScalar(1.0);
    xdot->PutScalar(2.0);
    SH->RecordSolution(1.0, *x, xdot.get());

    Teuchos::RCP<Epetra_Vector> y = Teuchos::rcp(new Epetra_Vector(*x));

    // interpolate 1st order
    SH->InterpolateSolution(0.5, *y, 2);

    double* norminf = new double[1];
    y->NormInf(norminf);
    CHECK_EQUAL(norminf[0], 0.25);

    // interpolate maximum order (3rd in this case)
    SH->InterpolateSolution(0.5, *y);

    y->NormInf(norminf);
    CHECK_EQUAL(norminf[0], 0.25);
  }
}
