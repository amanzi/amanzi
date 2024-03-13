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
#include "solver_fnbase6.hh"
#include "SolverNKA.hh"
#include "SolverNKA_LS.hh"
#include "SolverNKA_BT_ATS.hh"
#include "SolverNKA_LS_ATS.hh"
#include "SolverAA.hh"
#include "SolverLS.hh"
#include "SolverNox.hh"

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
  TEST_FIXTURE(test_data, NKA_SOLVER_EXACT_JACOBIAN_MONITOR_UPDATE)
  {
    std::cout << std::endl
              << "NKA solver, exact Jacobian (monitor update)..." << std::endl
              << "==============================================" << std::endl;

    // create the function class
    Teuchos::RCP<NonlinearProblem> fn = Teuchos::rcp(new NonlinearProblem(1.0, 1.0, true));

    // create the SolverState
    Teuchos::ParameterList plist;
    plist.set("nonlinear tolerance", 1e-8);
    plist.set("diverged tolerance", 1e10);
    plist.set("limit iterations", 10);
    plist.set("max du growth factor", 1e5);
    plist.set("max divergent iterations", 3);
    plist.set("max nka vectors", 1);
    plist.sublist("verbose object").set("verbosity level", "extreme");

    // create the Solver
    Teuchos::RCP<AmanziSolvers::SolverNKA<Epetra_Vector, Epetra_BlockMap>> nka =
      Teuchos::rcp(new AmanziSolvers::SolverNKA<Epetra_Vector, Epetra_BlockMap>(plist));
    nka->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.9;
    (*u)[1] = 0.9;

    // solve
    nka->Solve(u);
    CHECK_CLOSE(0.0, (*u)[0], 1.0e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.0e-6);

    // repeat solve with exact solution.
    std::cout << std::endl
              << "NKA solver, repeat solver..." << std::endl
              << "--------------------------------------------" << std::endl;
    nka->Solve(u);
    CHECK(nka->num_itrs() == 1);
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_SOLVER_EXACT_JACOBIAN_MONITOR_RESIDUAL)
  {
    std::cout << std::endl
              << "NKA solver, exact Jacobian (monitor residual)..." << std::endl
              << "==============================================" << std::endl;

    // create the function class
    Teuchos::RCP<NonlinearProblem> fn = Teuchos::rcp(new NonlinearProblem(1.0, 1.0, true));

    // create the SolverState
    Teuchos::ParameterList plist;
    plist.set("nonlinear tolerance", 1e-5);
    plist.set("diverged tolerance", 1e10);
    plist.set("limit iterations", 10);
    plist.set("max du growth factor", 1e5);
    plist.set("max divergent iterations", 3);
    plist.set("max nka vectors", 1);
    plist.set("monitor", "monitor residual");
    plist.set("make one iteration", true);
    plist.sublist("verbose object").set("verbosity level", "high");

    // create the Solver
    Teuchos::RCP<AmanziSolvers::SolverNKA<Epetra_Vector, Epetra_BlockMap>> nka =
      Teuchos::rcp(new AmanziSolvers::SolverNKA<Epetra_Vector, Epetra_BlockMap>(plist));
    nka->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.9;
    (*u)[1] = 0.9;

    // solve
    nka->Solve(u);
    CHECK_CLOSE(0.0, (*u)[0], 1.0e-5);
    CHECK_CLOSE(0.0, (*u)[1], 1.0e-5);

    // repeat solve with eaxt solution
    std::cout << std::endl
              << "NKA solver, repeat solver..." << std::endl
              << "--------------------------------------------" << std::endl;
    nka->Solve(u);
    CHECK(nka->num_itrs() == 1);
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_SOLVER_EXACT_JACOBIAN_GLOBALIZED)
  {
    std::cout << std::endl
              << "NKA solver, exact Jacobian (globalized)..." << std::endl
              << "============================================" << std::endl;

    // create the function class
    Teuchos::RCP<NonlinearProblem6> fn = Teuchos::rcp(new NonlinearProblem6(1.0, 1.0, true));

    // create the SolverState
    Teuchos::ParameterList plist;
    plist.set("nonlinear tolerance", 1e-8);
    plist.set("diverged tolerance", 1e10);
    plist.set("limit iterations", 10);
    plist.set("max du growth factor", 1e5);
    plist.set("max divergent iterations", 3);
    plist.set("max nka vectors", 1);
    plist.sublist("verbose object").set("verbosity level", "extreme");

    // create the Solver
    Teuchos::RCP<AmanziSolvers::SolverNKA<Epetra_Vector, Epetra_BlockMap>> nka =
      Teuchos::rcp(new AmanziSolvers::SolverNKA<Epetra_Vector, Epetra_BlockMap>(plist));
    nka->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.9;
    (*u)[1] = 0.9;

    // solve
    nka->Solve(u);
    CHECK_CLOSE(0.0, (*u)[0], 1.0e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.0e-6);
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_SOLVER_INEXACT_JACOBIAN)
  {
    std::cout << std::endl
              << "NKA solver, inexact Jacobian..." << std::endl
              << "==============================================" << std::endl;

    // create the function class
    Teuchos::RCP<NonlinearProblem> fn = Teuchos::rcp(new NonlinearProblem(1.0, 1.0, false));

    // create the SolverState
    Teuchos::ParameterList plist;
    plist.set("nonlinear tolerance", 1e-8);
    plist.set("diverged tolerance", 1e10);
    plist.set("limit iterations", 20);
    plist.set("max du growth factor", 1e5);
    plist.set("max divergent iterations", 3);
    plist.set("max nka vectors", 2);
    plist.sublist("verbose object").set("verbosity level", "high");

    // create the Solver
    Teuchos::RCP<AmanziSolvers::SolverNKA<Epetra_Vector, Epetra_BlockMap>> nka =
      Teuchos::rcp(new AmanziSolvers::SolverNKA<Epetra_Vector, Epetra_BlockMap>(plist));
    nka->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.9;
    (*u)[1] = 0.9;

    // solve
    nka->Solve(u);
    CHECK_CLOSE(0.0, (*u)[0], 1.e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.e-6);
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_SOLVER_INEXACT_JACOBIAN_GLOBALIZED)
  {
    std::cout << std::endl
              << "NKA solver, inexact Jacobian (globalized)..." << std::endl
              << "==============================================" << std::endl;


    // create the function class
    Teuchos::RCP<NonlinearProblem6> fn = Teuchos::rcp(new NonlinearProblem6(1.0, 1.0, false));

    // create the SolverState
    Teuchos::ParameterList plist;
    plist.set("nonlinear tolerance", 1e-8);
    plist.set("diverged tolerance", 1e10);
    plist.set("limit iterations", 20);
    plist.set("max du growth factor", 1e5);
    plist.set("max divergent iterations", 3);
    plist.set("max nka vectors", 2);
    plist.sublist("verbose object").set("verbosity level", "high");

    // create the Solver
    Teuchos::RCP<AmanziSolvers::SolverNKA<Epetra_Vector, Epetra_BlockMap>> nka =
      Teuchos::rcp(new AmanziSolvers::SolverNKA<Epetra_Vector, Epetra_BlockMap>(plist));
    nka->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.9;
    (*u)[1] = 0.9;

    // solve
    nka->Solve(u);
    CHECK_CLOSE(0.0, (*u)[0], 1.e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.e-6);
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_LS_SOLVER)
  {
    std::cout << std::endl
              << "NKA with line search..." << std::endl
              << "==============================================" << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverNKA_LS<Epetra_Vector, Epetra_BlockMap>> nka_bt =
      Teuchos::rcp(new AmanziSolvers::SolverNKA_LS<Epetra_Vector, Epetra_BlockMap>(plist));
    nka_bt->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.9;
    (*u)[1] = 0.9;

    // solve
    nka_bt->Solve(u);
    CHECK_CLOSE(0.0, (*u)[0], 1.0e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.0e-6);
  };

  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_LS_SOLVER_GLOBALIZATION)
  {
    std::cout << std::endl
              << "NKA with line search (globalized)..." << std::endl
              << "==============================================" << std::endl;

    // create the function class
    Teuchos::RCP<NonlinearProblem6> fn = Teuchos::rcp(new NonlinearProblem6(1.0, 1.0, true));

    // create the SolverState
    Teuchos::ParameterList plist;
    plist.set("nonlinear tolerance", 1e-6);
    plist.set("diverged tolerance", 1e10);
    plist.set("limit iterations", 15);
    plist.set("max du growth factor", 1e5);
    plist.set("max divergent iterations", 3);
    plist.sublist("verbose object").set("verbosity level", "high");

    // create the Solver
    Teuchos::RCP<AmanziSolvers::SolverNKA_LS<Epetra_Vector, Epetra_BlockMap>> nka_bt =
      Teuchos::rcp(new AmanziSolvers::SolverNKA_LS<Epetra_Vector, Epetra_BlockMap>(plist));
    nka_bt->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.9;
    (*u)[1] = 0.9;

    // solve
    nka_bt->Solve(u);
    CHECK_CLOSE(0.0, (*u)[0], 1.0e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.0e-6);
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_BT_ATS_SOLVER)
  {
    std::cout << std::endl
              << "NKA with backtracking, ATS custom..." << std::endl
              << "============================================" << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverNKA_BT_ATS<Epetra_Vector, Epetra_BlockMap>> nka_bt =
      Teuchos::rcp(new AmanziSolvers::SolverNKA_BT_ATS<Epetra_Vector, Epetra_BlockMap>(plist));
    nka_bt->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.9;
    (*u)[1] = 0.9;

    // solve
    nka_bt->Solve(u);

    std::cout << "Solution " << (*u)[0] << " " << (*u)[1] << "\n";

    CHECK_CLOSE(0.0, (*u)[0], 1.0e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.0e-6);
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_BT_ATS_SOLVER_GLOBALIZED)
  {
    std::cout << std::endl
              << "NKA with backtracking, ATS custom (globalized)..." << std::endl
              << "============================================" << std::endl;

    // create the function class
    Teuchos::RCP<NonlinearProblem6> fn = Teuchos::rcp(new NonlinearProblem6(1.0, 1.0, true));

    // create the SolverState
    Teuchos::ParameterList plist;
    plist.set("nonlinear tolerance", 1e-6);
    plist.set("diverged tolerance", 1e10);
    plist.set("limit iterations", 15);
    plist.set("max du growth factor", 1e5);
    plist.set("max divergent iterations", 3);
    plist.sublist("verbose object").set("verbosity level", "high");

    // create the Solver
    Teuchos::RCP<AmanziSolvers::SolverNKA_BT_ATS<Epetra_Vector, Epetra_BlockMap>> nka_bt =
      Teuchos::rcp(new AmanziSolvers::SolverNKA_BT_ATS<Epetra_Vector, Epetra_BlockMap>(plist));
    nka_bt->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.9;
    (*u)[1] = 0.9;

    // solve
    nka_bt->Solve(u);

    std::cout << "Solution " << (*u)[0] << " " << (*u)[1] << "\n";

    CHECK_CLOSE(0.0, (*u)[0], 1.0e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.0e-6);
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_LS_ATS_SOLVER)
  {
    std::cout << std::endl
              << "NKA with backtracking via line search, ATS custom..." << std::endl
              << "============================================" << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverNKA_LS_ATS<Epetra_Vector, Epetra_BlockMap>> nka_bt =
      Teuchos::rcp(new AmanziSolvers::SolverNKA_LS_ATS<Epetra_Vector, Epetra_BlockMap>(plist));
    nka_bt->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.9;
    (*u)[1] = 0.9;

    // solve
    nka_bt->Solve(u);

    std::cout << "Solution " << (*u)[0] << " " << (*u)[1] << "\n";

    CHECK_CLOSE(0.0, (*u)[0], 1.0e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.0e-6);
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_LS_ATS_SOLVER_GLOBALIZED)
  {
    std::cout << std::endl
              << "NKA with backtracking via line search, ATS custom (globalized)..." << std::endl
              << "============================================" << std::endl;

    // create the function class
    Teuchos::RCP<NonlinearProblem6> fn = Teuchos::rcp(new NonlinearProblem6(1.0, 1.0, true));

    // create the SolverState
    Teuchos::ParameterList plist;
    plist.set("nonlinear tolerance", 1e-6);
    plist.set("diverged tolerance", 1e10);
    plist.set("limit iterations", 15);
    plist.set("max du growth factor", 1e5);
    plist.set("max divergent iterations", 3);
    plist.sublist("verbose object").set("verbosity level", "high");

    // create the Solver
    Teuchos::RCP<AmanziSolvers::SolverNKA_LS_ATS<Epetra_Vector, Epetra_BlockMap>> nka_bt =
      Teuchos::rcp(new AmanziSolvers::SolverNKA_LS_ATS<Epetra_Vector, Epetra_BlockMap>(plist));
    nka_bt->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.9;
    (*u)[1] = 0.9;

    // solve
    nka_bt->Solve(u);

    std::cout << "Solution " << (*u)[0] << " " << (*u)[1] << "\n";

    CHECK_CLOSE(0.0, (*u)[0], 1.0e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.0e-6);
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, BT_LS_SOLVER)
  {
    std::cout << std::endl
              << "Line-search using Brent..." << std::endl
              << "============================================" << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverLS<Epetra_Vector, Epetra_BlockMap>> nka_bt =
      Teuchos::rcp(new AmanziSolvers::SolverLS<Epetra_Vector, Epetra_BlockMap>(plist));
    nka_bt->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.9;
    (*u)[1] = 0.9;

    // solve
    nka_bt->Solve(u);

    std::cout << "Solution " << (*u)[0] << " " << (*u)[1] << "\n";

    CHECK_CLOSE(0.0, (*u)[0], 1.0e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.0e-6);
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, BT_LS_SOLVER_GLOBALIZED)
  {
    std::cout << std::endl
              << "Line-search using Brent (globalized)..." << std::endl
              << "============================================" << std::endl;

    // create the function class
    Teuchos::RCP<NonlinearProblem6> fn = Teuchos::rcp(new NonlinearProblem6(1.0, 1.0, true));

    // create the SolverState
    Teuchos::ParameterList plist;
    plist.set("nonlinear tolerance", 1e-6);
    plist.set("diverged tolerance", 1e10);
    plist.set("limit iterations", 15);
    plist.set("max du growth factor", 1e5);
    plist.set("max divergent iterations", 3);
    plist.sublist("verbose object").set("verbosity level", "high");

    // create the Solver
    Teuchos::RCP<AmanziSolvers::SolverLS<Epetra_Vector, Epetra_BlockMap>> nka_bt =
      Teuchos::rcp(new AmanziSolvers::SolverLS<Epetra_Vector, Epetra_BlockMap>(plist));
    nka_bt->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.9;
    (*u)[1] = 0.9;

    // solve
    nka_bt->Solve(u);

    std::cout << "Solution " << (*u)[0] << " " << (*u)[1] << "\n";

    CHECK_CLOSE(0.0, (*u)[0], 1.0e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.0e-6);
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, AA_SOLVER)
  {
    std::cout << std::endl
              << "AA solver...." << std::endl
              << "============================================" << std::endl;

    // create the function class
    Teuchos::RCP<NonlinearProblem> fn = Teuchos::rcp(new NonlinearProblem(1.0, 1.0, true));

    // create the SolverState
    Teuchos::ParameterList plist;
    plist.set("nonlinear tolerance", 1e-7);
    plist.set("diverged tolerance", 1e10);
    plist.set("limit iterations", 15);
    plist.set("max du growth factor", 1e5);
    plist.set("max divergent iterations", 3);
    plist.set("max aa vectors", 4);
    plist.set("relaxation parameter", 1.);
    plist.sublist("verbose object").set("verbosity level", "high");

    // create the Solver
    Teuchos::RCP<AmanziSolvers::SolverAA<Epetra_Vector, Epetra_BlockMap>> aa =
      Teuchos::rcp(new AmanziSolvers::SolverAA<Epetra_Vector, Epetra_BlockMap>(plist));
    aa->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.95;
    (*u)[1] = 0.15;
    (*u)[2] = -0.51;
    (*u)[3] = 0.35;
    (*u)[4] = -0.54;

    // solve
    aa->Solve(u);
    CHECK_CLOSE(0.0, (*u)[0], 1.0e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.0e-6);
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, AA_SOLVER_GLOBALIZED)
  {
    std::cout << std::endl
              << "AA solver (globalized)...." << std::endl
              << "============================================" << std::endl;

    // create the function class
    Teuchos::RCP<NonlinearProblem6> fn = Teuchos::rcp(new NonlinearProblem6(1.0, 1.0, true));

    // create the SolverState
    Teuchos::ParameterList plist;
    plist.set("nonlinear tolerance", 1e-7);
    plist.set("diverged tolerance", 1e10);
    plist.set("limit iterations", 15);
    plist.set("max du growth factor", 1e5);
    plist.set("max divergent iterations", 3);
    plist.set("max aa vectors", 4);
    plist.set("relaxation parameter", 1.);
    plist.sublist("verbose object").set("verbosity level", "high");

    // create the Solver
    Teuchos::RCP<AmanziSolvers::SolverAA<Epetra_Vector, Epetra_BlockMap>> aa =
      Teuchos::rcp(new AmanziSolvers::SolverAA<Epetra_Vector, Epetra_BlockMap>(plist));
    aa->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.95;
    (*u)[1] = 0.15;
    (*u)[2] = -0.51;
    (*u)[3] = 0.35;
    (*u)[4] = -0.54;

    // solve
    aa->Solve(u);
    CHECK_CLOSE(0.0, (*u)[0], 1.0e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.0e-6);
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, NOX_SOLVER)
  {
    std::cout << std::endl
              << "NOX solver..." << std::endl
              << "============================================" << std::endl;

    // create the function class
    Teuchos::RCP<NonlinearProblem> fn = Teuchos::rcp(new NonlinearProblem(1.0, 1.0, true));

    // create the SolverState
    Teuchos::ParameterList plist;
    plist.sublist("nonlinear solver")
      .set("nonlinear tolerance", 1e-6)
      .set("diverged tolerance", 1e10)
      .set("limit iterations", 15)
      .set("max du growth factor", 1e5);
    plist.sublist("nonlinear solver").sublist("VerboseObject").set("Verbosity Level", "high");
    plist.sublist("JF matrix parameters");
    plist.sublist("linear operator").set("iterative method", "gmres");
    plist.sublist("linear operator").sublist("gmres parameters").set("size of Krylov space", 2);
    plist.sublist("linear operator").sublist("verbose object").set("verbosity level", "extreme");

    // create the Solver
    Teuchos::RCP<AmanziSolvers::SolverNox<Epetra_Vector, Epetra_BlockMap>> nox =
      Teuchos::rcp(new AmanziSolvers::SolverNox<Epetra_Vector, Epetra_BlockMap>(plist));
    nox->Init(fn, *map);

    // initial guess
    Teuchos::RCP<Epetra_Vector> u = Teuchos::rcp(new Epetra_Vector(*vec));
    (*u)[0] = -0.9;
    (*u)[1] = 0.9;

    // solve
    nox->Solve(u);
    CHECK_CLOSE(0.0, (*u)[0], 1.0e-6);
    CHECK_CLOSE(0.0, (*u)[1], 1.0e-6);
  };
} // SUITE
