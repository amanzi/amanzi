#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "solver_fnbase1_cont.hh"
#include "SolverContinuation.hh"

using namespace Amanzi;

SUITE(SOLVERS) {
// data structures for testing
struct test_data {
  Teuchos::RCP<Map_type> map;
  Teuchos::RCP<Vector_type> vec;

  test_data() {
    auto comm = getDefaultComm(); 
    map = Teuchos::rcp(new Epetra_Map(2, 0, *comm));
    vec = Teuchos::rcp(new Epetra_Vector(*map));
  }

};


/* ******************************************************************/
TEST_FIXTURE(test_data, CONT_SOLVER_EXACT_JACOBIAN) {
  std::cout << "Continuation nonlinear solver, exact Jacobian..." << std::endl;

  // create the function class
  Teuchos::RCP<NonlinearProblem> fn = Teuchos::rcp(new NonlinearProblem(1.0, 1.0, true));

  // create the SolverState
  Teuchos::ParameterList list;
  Teuchos::ParameterList& inner = list.sublist("inner solver");
  inner.set("solver type", "nka");
  Teuchos::ParameterList& plist = inner.sublist("nka parameters");
  plist.set("nonlinear tolerance", 1e-8);
  plist.set("diverged tolerance", 1e10);
  plist.set("limit iterations", 10);
  plist.set("max du growth factor", 1e5);
  plist.set("max divergent iterations", 3);
  plist.set("max nka vectors", 1);
  plist.sublist("verbose object").set("verbosity level", "extreme");

  // create the Solver
  auto cont =
      Teuchos::rcp(new AmanziSolvers::SolverContinuation<Epetra_Vector, Epetra_BlockMap>(list));
  cont->Init(fn, *map);

  // initial guess
  Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
  (*u)[0] = -0.9;
  (*u)[1] =  0.9; 

  // solve
  cont->Solve(u);
  CHECK_CLOSE(0.7937005259840998, (*u)[0], 1.0e-6);
  CHECK_CLOSE(0.7937005259840998, (*u)[1], 1.0e-6);
};


/* ******************************************************************/
TEST_FIXTURE(test_data, CONT_SOLVER_INEXACT_JACOBIAN) {
  std::cout << std::endl 
            << "Continuation nonlinear solver, inexact Jacobian..." << std::endl;

  // create the function class
  Teuchos::RCP<NonlinearProblem> fn = Teuchos::rcp(new NonlinearProblem(1.0, 1.0, false));

  // create the SolverState
  Teuchos::ParameterList list;
  Teuchos::ParameterList& inner = list.sublist("inner solver");
  inner.set("solver type", "nka");
  Teuchos::ParameterList& plist = inner.sublist("nka parameters");
  plist.set("nonlinear tolerance", 1e-8);
  plist.set("diverged tolerance", 1e10);
  plist.set("limit iterations", 20);
  plist.set("max du growth factor", 1e5);
  plist.set("max divergent iterations", 3);
  plist.set("max nka vectors", 2);
  plist.sublist("verbose object").set("verbosity level", "high");

  // create the Solver
  Teuchos::RCP<AmanziSolvers::SolverContinuation<Epetra_Vector, Epetra_BlockMap> > cont =
      Teuchos::rcp(new AmanziSolvers::SolverContinuation<Epetra_Vector, Epetra_BlockMap>(list));
  cont->Init(fn, *map);

  // initial guess
  Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
  (*u)[0] = -0.9;
  (*u)[1] =  0.9;

  // solve
  cont->Solve(u);
  CHECK_CLOSE(0.7937005259840998, (*u)[0], 1.0e-6);
  CHECK_CLOSE(0.7937005259840998, (*u)[1], 1.0e-6);
};
}  // SUITE




