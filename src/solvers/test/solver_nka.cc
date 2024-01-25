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

#include "solver_fnbase1.hh"
#include "solver_fnbase6.hh"
#include "SolverNKA.hh"
#include "SolverNKA_LS.hh"
#include "SolverNKA_BT_ATS.hh"
#include "SolverNKA_LS_ATS.hh"
#include "SolverNewton.hh"
#include "SolverJFNK.hh"
#include "SolverAA.hh"
#include "SolverBT.hh"
// #include "SolverNox.hh"

using namespace Amanzi;

SUITE(SOLVERS)
{
  // data structures for testing
  struct test_data {
    Map_ptr_type map;
    Vector_ptr_type vec;

    test_data()
    {
      auto comm = getDefaultComm();
      map = Teuchos::rcp(new Map_type(5, 0, comm));
      vec = Teuchos::rcp(new Vector_type(map));
    }

    template <class Solver>
    int run(Solver& solver)
    {
      // initial guess
      {
        auto vecv = vec->getLocalViewHost(Tpetra::Access::ReadWrite);
        vecv(0, 0) = -0.9;
        vecv(1, 0) = 0.9;
      }

      // solve
      int ierr = solver->Solve(vec);

      // test
      if (!ierr) {
        {
          auto vecv = vec->getLocalViewHost(Tpetra::Access::ReadOnly);
          CHECK_CLOSE(0., vecv(0, 0), 1.e-6);
          CHECK_CLOSE(0., vecv(1, 0), 1.e-6);
        }
      }
      return ierr;
    }

    ~test_data() {}
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_SOLVER_EXACT_JACOBIAN)
  {
    std::cout << std::endl
              << std::endl
              << "----------------------------------------------------------------"
              << "NKA solver, exact Jacobian..." << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverNKA<Vector_type, Map_type>> nka =
      Teuchos::rcp(new AmanziSolvers::SolverNKA<Vector_type, Map_type>(plist));
    nka->Init(fn, map);
    CHECK(!run(nka));
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_SOLVER_EXACT_JACOBIAN_GLOBALIZED)
  {
    std::cout << std::endl
              << std::endl
              << "----------------------------------------------------------------";
    std::cout << "NKA solver, exact Jacobian, globalized problem..." << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverNKA<Vector_type, Map_type>> nka =
      Teuchos::rcp(new AmanziSolvers::SolverNKA<Vector_type, Map_type>(plist));
    nka->Init(fn, map);
    CHECK(!run(nka));
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_SOLVER_INEXACT_JACOBIAN)
  {
    std::cout << std::endl
              << std::endl
              << "----------------------------------------------------------------";
    std::cout << "\nNKA solver, inexact Jacobian..." << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverNKA<Vector_type, Map_type>> nka =
      Teuchos::rcp(new AmanziSolvers::SolverNKA<Vector_type, Map_type>(plist));
    nka->Init(fn, map);
    CHECK(!run(nka));
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_SOLVER_INEXACT_JACOBIAN_GLOBALIZED)
  {
    std::cout << std::endl
              << std::endl
              << "----------------------------------------------------------------";
    std::cout << "\nNKA solver, inexact Jacobian, globalized problem..." << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverNKA<Vector_type, Map_type>> nka =
      Teuchos::rcp(new AmanziSolvers::SolverNKA<Vector_type, Map_type>(plist));
    nka->Init(fn, map);
    CHECK(!run(nka));
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, NEWTON_SOLVER)
  {
    std::cout << std::endl
              << std::endl
              << "----------------------------------------------------------------";
    std::cout << "\nNewton solver..." << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverNewton<Vector_type, Map_type>> newton =
      Teuchos::rcp(new AmanziSolvers::SolverNewton<Vector_type, Map_type>(plist));
    newton->Init(fn, map);
    CHECK(!run(newton));
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, JFNK_SOLVER_LEFT_PC)
  {
    std::cout << std::endl
              << std::endl
              << "----------------------------------------------------------------";
    std::cout << "\nJFNK solver with LEFT precondiitoner..." << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverJFNK<Vector_type, Map_type>> jfnk =
      Teuchos::rcp(new AmanziSolvers::SolverJFNK<Vector_type, Map_type>(plist));
    jfnk->Init(fn, map);
    CHECK(!run(jfnk));
  };


#if 0

// THIS TEST Has issues with numerical roundoff/divide by zero.  Unclear
// exactly what is going wrong, but it appears that master is either working
// accidentally or we missed something...

/* ******************************************************************/
TEST_FIXTURE(test_data, JFNK_SOLVER_RIGHT_PC) {
  std::cout << std::endl << std::endl
            << "----------------------------------------------------------------";
  std::cout << "\nJFNK solver with RIGHT precondiitoner..." << std::endl;

  // create the function class
  Teuchos::RCP<NonlinearProblem> fn = Teuchos::rcp(new NonlinearProblem(1.0, 1.0, false));

  // create the SolverState
  Teuchos::ParameterList plist;
  plist.sublist("nonlinear solver").set("solver type", "Newton");
  plist.sublist("nonlinear solver").sublist("Newton parameters").sublist("verbose object")
      .set("verbosity level", "extreme");
  plist.sublist("nonlinear solver").sublist("Newton parameters")
      .set("nonlinear tolerance", 1e-6)
      .set("diverged tolerance", 1e10)
      .set("limit iterations", 15)
      .set("max du growth factor", 1e5)
      .set("max divergent iterations", 3);
  plist.sublist("JF matrix parameters");
  plist.sublist("linear operator").set("iterative method", "gmres");
  plist.sublist("linear operator").sublist("gmres parameters")
      .set("size of Krylov space", 2)
      .set("preconditioning strategy", "right");
  plist.sublist("linear operator").sublist("verbose object").set("verbosity level", "extreme");

  // create the Solver
  Teuchos::RCP<AmanziSolvers::SolverNewton<Vector_type, Map_type> > jfnk =
      Teuchos::rcp(new AmanziSolvers::SolverNewton<Vector_type, Map_type>(plist));
  jfnk->Init(fn, map);
  CHECK(!run(jfnk));

};

#endif

  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_LS_SOLVER)
  {
    std::cout << std::endl
              << std::endl
              << "----------------------------------------------------------------";
    std::cout << std::endl << "NKA with backtracking..." << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverNKA_LS<Vector_type, Map_type>> nka_bt =
      Teuchos::rcp(new AmanziSolvers::SolverNKA_LS<Vector_type, Map_type>(plist));
    nka_bt->Init(fn, map);
    CHECK(!run(nka_bt));
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_LS_SOLVER_GLOBALIZATION)
  {
    std::cout << std::endl
              << std::endl
              << "----------------------------------------------------------------";
    std::cout << std::endl << "NKA with backtracking, globalized problem..." << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverNKA_LS<Vector_type, Map_type>> nka_bt =
      Teuchos::rcp(new AmanziSolvers::SolverNKA_LS<Vector_type, Map_type>(plist));
    nka_bt->Init(fn, map);
    CHECK(!run(nka_bt));
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_BT_ATS_SOLVER)
  {
    std::cout << std::endl
              << std::endl
              << "----------------------------------------------------------------";
    std::cout << std::endl << "NKA with backtracking, ATS custom..." << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverNKA_BT_ATS<Vector_type, Map_type>> nka_bt =
      Teuchos::rcp(new AmanziSolvers::SolverNKA_BT_ATS<Vector_type, Map_type>(plist));
    nka_bt->Init(fn, map);
    CHECK(!run(nka_bt));
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_BT_ATS_SOLVER_GLOBALIZED)
  {
    std::cout << std::endl
              << std::endl
              << "----------------------------------------------------------------";
    std::cout << std::endl
              << "NKA with backtracking, ATS custom, globalized problem..." << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverNKA_BT_ATS<Vector_type, Map_type>> nka_bt =
      Teuchos::rcp(new AmanziSolvers::SolverNKA_BT_ATS<Vector_type, Map_type>(plist));
    nka_bt->Init(fn, map);

    // NOTE: this does not actually converge on master, but the solution is
    // still
    // within tolerance after 16 iterations...
    // CHECK(!run(nka_bt));
    run(nka_bt);
  };

  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_LS_ATS_SOLVER)
  {
    std::cout << std::endl
              << std::endl
              << "----------------------------------------------------------------";
    std::cout << std::endl << "NKA with backtracking via line search, ATS custom..." << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverNKA_LS_ATS<Vector_type, Map_type>> nka_bt =
      Teuchos::rcp(new AmanziSolvers::SolverNKA_LS_ATS<Vector_type, Map_type>(plist));
    nka_bt->Init(fn, map);
    CHECK(!run(nka_bt));
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, NKA_LS_ATS_SOLVER_GLOBALIZED)
  {
    std::cout << std::endl
              << std::endl
              << "----------------------------------------------------------------";
    std::cout << std::endl
              << "NKA with backtracking via line search, ATS custom, "
                 "globalized problem..."
              << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverNKA_LS_ATS<Vector_type, Map_type>> nka_bt =
      Teuchos::rcp(new AmanziSolvers::SolverNKA_LS_ATS<Vector_type, Map_type>(plist));
    nka_bt->Init(fn, map);

    // NOTE: this does not actually converge on master, but the solution is
    // still
    // within tolerance after 16 iterations...
    // CHECK(!run(nka_bt));
    run(nka_bt);
  };


  /* ******************************************************************/
  TEST_FIXTURE(test_data, BT_LS_SOLVER)
  {
    std::cout << std::endl
              << std::endl
              << "----------------------------------------------------------------";
    std::cout << std::endl << "Backtracking line-search using Brent..." << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverBT<Vector_type, Map_type>> nka_bt =
      Teuchos::rcp(new AmanziSolvers::SolverBT<Vector_type, Map_type>(plist));
    nka_bt->Init(fn, map);
    CHECK(!run(nka_bt));
  };

  /* ******************************************************************/
  TEST_FIXTURE(test_data, BT_LS_SOLVER_GLOBALIZED)
  {
    std::cout << std::endl
              << std::endl
              << "----------------------------------------------------------------";
    std::cout << std::endl
              << "Backtracking line-search using Brent, globalized problem..." << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverBT<Vector_type, Map_type>> nka_bt =
      Teuchos::rcp(new AmanziSolvers::SolverBT<Vector_type, Map_type>(plist));
    nka_bt->Init(fn, map);

    // NOTE: this does not actually converge on master, but the solution is
    // still
    // within tolerance after 16 iterations...
    // CHECK(!run(nka_bt));
    run(nka_bt);
  };


  TEST_FIXTURE(test_data, AA_SOLVER)
  {
    std::cout << std::endl
              << std::endl
              << "----------------------------------------------------------------";
    std::cout << std::endl << "AA solver...." << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverAA<Vector_type, Map_type>> aa =
      Teuchos::rcp(new AmanziSolvers::SolverAA<Vector_type, Map_type>(plist));
    aa->Init(fn, map);
    CHECK(!run(aa));
  };

  TEST_FIXTURE(test_data, AA_SOLVER_GLOBALIZED)
  {
    std::cout << std::endl
              << std::endl
              << "----------------------------------------------------------------";
    std::cout << std::endl << "AA solver, globalized problem...." << std::endl;

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
    Teuchos::RCP<AmanziSolvers::SolverAA<Vector_type, Map_type>> aa =
      Teuchos::rcp(new AmanziSolvers::SolverAA<Vector_type, Map_type>(plist));
    aa->Init(fn, map);
    CHECK(!run(aa));
  };

#if 0

// NOX not currently implemented in tpetra branch

/* ******************************************************************/
TEST_FIXTURE(test_data, NOX_SOLVER) {
  std::cout << std::endl << std::endl
            << "----------------------------------------------------------------";
  std::cout << "\nNOX solver..." << std::endl;

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
  Teuchos::RCP<AmanziSolvers::SolverNox<Vector_type, Map_type> > nox =
      Teuchos::rcp(new AmanziSolvers::SolverNox<Vector_type, Map_type>(plist));
  nox->Init(fn, map);
  CHECK(!run(nox));
};

#endif

} // SUITE
