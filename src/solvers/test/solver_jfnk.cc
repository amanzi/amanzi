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
#include "Epetra_Vector.h"

#include "solver_fnbase1.hh"
#include "SolverJFNK.hh"

using namespace Amanzi;

SUITE(SOLVERS)
{
  // data structures for testing
  struct test_data {
    Epetra_MpiComm* comm;
    Teuchos::RCP<Epetra_Map> map;
    Teuchos::RCP<Epetra_Vector> vec;

    test_data()
    {
      comm = new Epetra_MpiComm(MPI_COMM_SELF);
      map = Teuchos::rcp(new Epetra_Map(5, 0, *comm));
      vec = Teuchos::rcp(new Epetra_Vector(*map));
    }

    ~test_data() { delete comm; }
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, JFNK_SOLVER_LEFT_PC)
  {
    std::cout << "JFNK solver with LEFT precondiitoner..." << std::endl;

    // create the function class
    Teuchos::RCP<NonlinearProblem> fn = Teuchos::rcp(new NonlinearProblem(1.0, 1.0, false));

    // create the SolverState
    Teuchos::ParameterList plist;
    plist.sublist("nonlinear solver").set("solver type", "Newton");
    plist.sublist("nonlinear solver")
      .sublist("Newton parameters")
      .sublist("verbose object")
      .set("verbosity level", "extreme");
    plist.sublist("nonlinear solver")
      .sublist("Newton parameters")
      .set("nonlinear tolerance", 1e-6)
      .set("diverged tolerance", 1e10)
      .set("limit iterations", 15)
      .set("max du growth factor", 1e5)
      .set("max divergent iterations", 3);
    plist.sublist("JF matrix parameters");
    plist.sublist("linear operator").set("iterative method", "gmres");
    plist.sublist("linear operator").sublist("gmres parameters").set("size of Krylov space", 2);
    plist.sublist("linear operator").sublist("verbose object").set("verbosity level", "extreme");

    // create the Solver
    Teuchos::RCP<AmanziSolvers::SolverJFNK<Epetra_Vector, Epetra_BlockMap>> jfnk =
      Teuchos::rcp(new AmanziSolvers::SolverJFNK<Epetra_Vector, Epetra_BlockMap>(plist));
    jfnk->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.9;
    (*u)[1] = 0.9;

    // solve
    jfnk->Solve(u);
    CHECK_CLOSE(0.0, (*u)[0], 1.0e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.0e-6);
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, JFNK_SOLVER_RIGHT_PC)
  {
    std::cout << "\nJFNK solver with RIGHT preconditioner..." << std::endl;

    // create the function class
    Teuchos::RCP<NonlinearProblem> fn = Teuchos::rcp(new NonlinearProblem(1.0, 1.0, false));

    // create the SolverState
    Teuchos::ParameterList plist;
    Teuchos::Array<std::string> criteria;
    criteria.push_back("relative rhs");
    criteria.push_back("absolute residual");

    plist.sublist("nonlinear solver").set("solver type", "Newton");
    plist.sublist("nonlinear solver")
      .sublist("Newton parameters")
      .sublist("verbose object")
      .set("verbosity level", "extreme");
    plist.sublist("nonlinear solver")
      .sublist("Newton parameters")
      .set("nonlinear tolerance", 1e-6)
      .set("diverged tolerance", 1e10)
      .set("limit iterations", 15)
      .set("max du growth factor", 1e5)
      .set("max divergent iterations", 3);
    plist.sublist("JF matrix parameters");
    plist.sublist("linear operator").set("iterative method", "gmres");
    plist.sublist("linear operator")
      .sublist("gmres parameters")
      .set("size of Krylov space", 2)
      .set("preconditioning strategy", "right")
      .set<double>("error tolerance", 1e-14)
      .set<Teuchos::Array<std::string>>("convergence criteria", criteria);
    plist.sublist("linear operator").sublist("verbose object").set("verbosity level", "extreme");

    // create the Solver
    Teuchos::RCP<AmanziSolvers::SolverJFNK<Epetra_Vector, Epetra_BlockMap>> jfnk =
      Teuchos::rcp(new AmanziSolvers::SolverJFNK<Epetra_Vector, Epetra_BlockMap>(plist));
    jfnk->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.9;
    (*u)[1] = 0.9;

    // solve
    jfnk->Solve(u);
    CHECK_CLOSE(0.0, (*u)[0], 1.0e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.0e-6);
  };
} // SUITE
