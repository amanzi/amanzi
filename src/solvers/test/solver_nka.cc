#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"

#include "solver_fnbase6.hh"
#include "SolverNKA.hh"

using namespace Amanzi;

SUITE(SOLVERS) {
  // data structures for testing
  struct test_data {
    Epetra_MpiComm *comm;
    Teuchos::RCP<Epetra_Map> map;
    Teuchos::RCP<Epetra_Vector> vec;

    test_data() {
      comm = new Epetra_MpiComm(MPI_COMM_SELF);
      map = Teuchos::rcp(new Epetra_Map(2,0,*comm));
      vec = Teuchos::rcp(new Epetra_Vector(*map));
    }

    ~test_data() { delete comm; }
  };

  TEST_FIXTURE(test_data, NKA_SOLVER_EXACT_JACOBIAN) {
    std::cout << "NKA nonlinear solver, exact Jacobian..." << std::endl;

    // create the function class
    Teuchos::RCP<NonlinearProblem> fn = Teuchos::rcp(new NonlinearProblem(1.0, 1.0, true));

    // create the SolverState
    Teuchos::ParameterList plist;
    plist.set("nonlinear tolerance", 1.e-6);
    plist.set("diverged tolerance", 1.e10);
    plist.set("limit iterations", 100);
    plist.set("max du growth factor", 1.0e5);
    plist.set("max divergent iterations", 3);
    plist.set("max nka vectors", 1);
    plist.set("convergence monitor", "monitor update");
    plist.sublist("VerboseObject").set("Verbosity Level", "extreme");

    // create the Solver
    Teuchos::RCP<Amanzi::AmanziSolvers::SolverNKA<Epetra_Vector> > nka =
        Teuchos::rcp(new Amanzi::AmanziSolvers::SolverNKA<Epetra_Vector>(plist, fn, *vec));

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.9;
    (*u)[1] =  0.9;

    // solve
    nka->Solve(u);
    CHECK_CLOSE(0.0, (*u)[0], 1.0e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.0e-6);
  };


  TEST_FIXTURE(test_data, NKA_SOLVER_INEXACT_JACOBIAN) {
    std::cout << "NKA nonlinear solver, inexact Jacobian..." << std::endl;

    // create the function class
    Teuchos::RCP<NonlinearProblem> fn = Teuchos::rcp(new NonlinearProblem(1.0, 1.0, false));

    // create the SolverState
    Teuchos::ParameterList plist;
    plist.set("nonlinear tolerance", 1.e-6);
    plist.set("diverged tolerance", 1.e10);
    plist.set("limit iterations", 100);
    plist.set("max du growth factor", 1.e5);
    plist.set("max divergent iterations", 3);
    plist.set("max nka vectors", 3);
    plist.set("convergence monitor", "monitor update");
    plist.sublist("VerboseObject").set("Verbosity Level", "extreme");

    // create the Solver
    Teuchos::RCP<Amanzi::AmanziSolvers::SolverNKA<Epetra_Vector> > nka =
        Teuchos::rcp(new Amanzi::AmanziSolvers::SolverNKA<Epetra_Vector>(plist, fn, *vec));

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.9;
    (*u)[1] =  0.9;

    // solve
    nka->Solve(u);
    CHECK_CLOSE(0.0, (*u)[0], 1.e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.e-6);
  };
}




