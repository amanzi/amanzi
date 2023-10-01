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
#include "SolverNewton.hh"

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
  TEST_FIXTURE(test_data, NEWTON_SOLVER)
  {
    std::cout << "Newton solver..." << std::endl;

    // create the function class
    Teuchos::RCP<NonlinearProblem> fn = Teuchos::rcp(new NonlinearProblem(1.0, 1.0, true));

    // create the SolverState
    Teuchos::ParameterList plist;
    plist.set("nonlinear tolerance", 1e-6);
    plist.set("diverged tolerance", 1e10);
    plist.set("limit iterations", 15);
    plist.set("max du growth factor", 1e5);
    plist.set("max divergent iterations", 3);
    plist.sublist("verbose object").set("verbosity level", "high");

    // create the Solver
    auto newton = Teuchos::rcp(new AmanziSolvers::SolverNewton<Epetra_Vector, Epetra_BlockMap>(plist));
    newton->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.9;
    (*u)[1] = 0.9;

    // solve
    newton->Solve(u);
    CHECK_CLOSE(0.0, (*u)[0], 1.0e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.0e-6);

    // verify quadratic convergence
    double r;
    const auto& history = newton->history();
    int n = history.size();
    for (int i = 2; i < n; ++i) {
      r = std::log(history[i].first) / std::log(history[i - 1].first);
      CHECK(r > 2.0);
      r = std::log(history[i].second) / std::log(history[i - 1].second);
      CHECK(r > 2.0);
    }

    // repeat solve with exact solution. 
    std::cout << "\nNetown solver, repeat solver..." << std::endl;
    newton->Solve(u);
    CHECK(newton->num_itrs() == 1);
  };
} // SUITE
