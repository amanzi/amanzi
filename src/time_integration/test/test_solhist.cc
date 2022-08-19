#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Vector.h"

#include "SolutionHistory.hh"

using namespace Amanzi;

SUITE(SolutionHistoryTests) {

  // data structures for testing
  struct test_data {
    Epetra_MpiComm *comm;
    Teuchos::RCP<Epetra_Vector> x;

    test_data() {
      comm = new Epetra_MpiComm(MPI_COMM_SELF);
      Epetra_Map map(2,0,*comm);
      x = Teuchos::rcp(new Epetra_Vector(map));
    }
    ~test_data() {
      delete comm;
    }
  };

  TEST_FIXTURE(test_data, SolutionHistory_1) {
    std::cout << "Test: SolutionHistory_1" << std::endl;

    // create a solution history of size three
    Teuchos::RCP<Amanzi::SolutionHistory<Epetra_Vector> > SH =
      Teuchos::rcp(new Amanzi::SolutionHistory<Epetra_Vector>("myhist", 3, 0.0, *x));
    x->PutScalar(1.0);
    CHECK_EQUAL(SH->history_size(), 1);
    CHECK_EQUAL(SH->MostRecentTime(), 0.0);

    SH->RecordSolution(1.0,* x);

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
    Teuchos::RCP<Epetra_Vector> y = Teuchos::rcp(new Epetra_Vector(*x));
    SH->MostRecentSolution(*y);

    double *norminf = new double[1];
    y->NormInf(norminf); // compute the max norm of the computed difference
    CHECK_EQUAL(*norminf, 9.0);

    std::vector<double> h;
    SH->TimeDeltas(h);

    CHECK_EQUAL(h.size(), SH->history_size() - 1);
    CHECK_EQUAL(h[0], 1.0);
    CHECK_EQUAL(h[1], 2.0);

    // interpolate 1st order
    SH->InterpolateSolution(2.5, *y, 1);

    // we should get 6.5
    y->NormInf(norminf);
    CHECK_EQUAL(norminf[0], 6.5);

    // interpolate 2nd order
    SH->InterpolateSolution(2.5, *y, 2);

    // from the quadratic through (1,1),(2,4),(3,9)
    // we should get 6.25
    // x->PutScalar(6.25);
    // x->Update(-1.0, *y, 1.0);
    y->NormInf(norminf);
    CHECK_EQUAL(norminf[0], 6.25);

    // interpolate maximum order (should be 2nd)
    SH->InterpolateSolution(2.5, *y);

    // from the quadratic through (1,1),(2,4),(3,9)
    // we should get 6.25
    y->NormInf(norminf);
    CHECK_EQUAL(norminf[0], 6.25);
  }

  TEST_FIXTURE(test_data, SolutionHistory_2) {
    std::cout << "Test: SolutionHistory_2" << std::endl;
    Teuchos::RCP<Epetra_Vector> xdot = Teuchos::rcp(new Epetra_Vector(*x));

    x->PutScalar(0.0);
    xdot->PutScalar(0.0);
    // create a solution history of size three
    Teuchos::RCP<Amanzi::SolutionHistory<Epetra_Vector> > SH =
      Teuchos::rcp(new Amanzi::SolutionHistory<Epetra_Vector>("myhist", 4, 0.0, *x, xdot.get()));

    x->PutScalar(1.0);
    xdot->PutScalar(2.0);
    SH->RecordSolution(1.0, *x, xdot.get());

    Teuchos::RCP<Epetra_Vector> y = Teuchos::rcp(new Epetra_Vector(*x));

    // interpolate 1st order
    SH->InterpolateSolution(0.5, *y, 2);

    double *norminf = new double[1];
    y->NormInf(norminf);
    CHECK_EQUAL(norminf[0], 0.25);

    // interpolate maximum order (3rd in this case)
    SH->InterpolateSolution(0.5, *y);

    y->NormInf(norminf);
    CHECK_EQUAL(norminf[0], 0.25);
  }
}
