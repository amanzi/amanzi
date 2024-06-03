/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//!
#include <ios>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerbosityLevel.hpp"

#include "dbc.hh"
#include "exceptions.hh"
#include "errors.hh"

#include "AmanziTypes.hh"
#include "AmanziComm.hh"
#include "AmanziVector.hh"

#include "BDFFnBase.hh"
#include "BDF1_TI.hh"
#include "FnBaseDefs.hh"

#include "ode_fnbase.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziSolvers;


SUITE(ODEIntegrationTests)
{
  // data structures for testing
  struct test_data {
    Comm_ptr_type comm;
    Vector_ptr_type init, u, u_dot, u_ex;
    Teuchos::RCP<Teuchos::ParameterList> plist;

    test_data()
    {
      // comm and mesh for maps
      comm = getDefaultComm();
      plist = Teuchos::rcp(new Teuchos::ParameterList());

      auto map = Teuchos::rcp(new Map_type(2, 0, comm));

      init = Teuchos::rcp(new Vector_type(map));
      u = Teuchos::rcp(new Vector_type(map));
      u_dot = Teuchos::rcp(new Vector_type(map));
      u_ex = Teuchos::rcp(new Vector_type(map));
    }
  };


  TEST_FIXTURE(test_data, NonlinearODE_BDF1_NKA_TrueJacobian)
  {
    std::cout << "Test: NonlinearODE_bdf1 with NKA, PC is True Jacobian" << std::endl;

    // set the parameter list for BDF1
    // set the parameter list for BDF1
    // set up solver params
    plist->set("solver type", "nka");
    plist->sublist("nka parameters").set("limit iterations", 20);
    plist->sublist("nka parameters").set("nonlinear tolerance", 1e-10);
    plist->sublist("nka parameters").set("diverged tolerance", 1.0e4);
    plist->sublist("nka parameters").set("convergence monitor", "monitor update");

    // set up time integrator params
    plist->set("timestep controller type", "standard");
    plist->sublist("timestep controller standard parameters")
      .set("time step reduction factor", 0.9);
    plist->sublist("timestep controller standard parameters").set("time step increase factor", 1.1);
    plist->sublist("timestep controller standard parameters").set("max time step", 5e-3);
    plist->sublist("timestep controller standard parameters").set("min time step", 1e-20);
    plist->sublist("timestep controller standard parameters").set("min iterations", 5);
    plist->sublist("timestep controller standard parameters").set("max iterations", 10);
    plist->sublist("timestep controller standard parameters")
      .set("preconditioner lag iterations", 2);

    // create the PDE problem
    nonlinearODE NF(1., 1., true);

    // create the time stepper
    Teuchos::RCP<Amanzi::BDF1_TI<Vector_type, Map_type>> TS =
      Teuchos::rcp(new BDF1_TI<Vector_type, Map_type>(NF, plist, init));

    // initial value
    u->putScalar(-1.0);
    u_dot->putScalar(1.0);

    // initial time
    double t = 0.0;

    // final time
    double tout = 2.0;

    // initial time step
    double h = 1.0e-5;
    double hnext;
    // initialize the state of the time stepper
    TS->SetInitialState(t, u, u_dot);

    // iterate until the final time
    int i = 0;
    double tlast = t;

    std::cout << "starting time integration" << std::endl;
    do {
      if (tlast + h > tout) {
        std::cout << "adjusting h, to hit the final time exactly...\n";
        h = tout - tlast;
      }

      bool redo(false);
      do {
        redo = TS->TimeStep(h, hnext, u);
      } while (redo);

      u->print(std::cout);

      TS->CommitSolution(h, u);

      h = hnext;
      i++;

      tlast = TS->time();
    } while (tout > tlast);

    // compute the error with the exact solution
    u_ex->putScalar(-1.0 / 3.0);
    u->update(1.0, *u_ex, -1.0);

    double norm;
    norm = u->normInf();

    TS->ReportStatistics_();
    CHECK_CLOSE(0.0, norm, 1e-3);
  }

  // TEST_FIXTURE(test_data, NonlinearODE_BDF1_NKA_ApproxJacobian) {
  //   std::cout << "Test: NonlinearODE_bdf1 with NKA, PC is 1/h Jacobian" <<
  //   std::endl;

  //   // set the parameter list for BDF1
  //   // set up solver params
  //   plist->set("solver type", "nka");
  //   plist->sublist("nka parameters").set("limit iterations", 20);
  //   plist->sublist("nka parameters").set("nonlinear tolerance",1e-10);
  //   plist->sublist("nka parameters").set("diverged tolerance",1.0e4);
  //   plist->sublist("nka parameters").set("convergence monitor","monitor
  //   update");

  //   // set up time integrator params
  //   plist->set("timestep controller type", "standard");
  //   plist->sublist("timestep controller standard parameters").set("time step
  //   reduction factor",0.9); plist->sublist("timestep controller standard
  //   parameters").set("time step increase factor",1.1);
  //   plist->sublist("timestep controller standard parameters").set("max time
  //   step",5e-3); plist->sublist("timestep controller standard
  //   parameters").set("min time step",1e-20); plist->sublist("timestep
  //   controller standard parameters").set("min iterations", 5);
  //   plist->sublist("timestep controller standard parameters").set("max
  //   iterations",10); plist->sublist("timestep controller standard
  //   parameters").set("preconditioner lag iterations",2);

  //   // create the PDE problem
  //   nonlinearODE NF(1., 1., false);

  //   // create the time stepper
  //   Teuchos::RCP<Amanzi::BDF1_TI<Vector_type, Epetra_BlockMap> > TS =
  //       Teuchos::rcp(new BDF1_TI<Vector_type, Epetra_BlockMap>(NF, plist,
  //       init));

  //   // initial value
  //   u->putScalar(-1.0);
  //   u_dot->putScalar(1.0);

  //   // initial time
  //   double t=0.0;

  //   // final time
  //   double tout = 2.0;

  //   // initial time step
  //   double h = 1.0e-5;
  //   double hnext;
  //   // initialize the state of the time stepper
  //   TS->SetInitialState(t, u, u_dot);

  //   // iterate until the final time
  //   int i=0;
  //   double tlast = t;

  //   std::cout << "starting time integration" << std::endl;
  //   do {
  //     if (tlast + h > tout) {
  //       std::cout << "adjusting h, to hit the final time exactly...\n";
  //       h = tout - tlast;
  //     }

  //     bool redo(false);
  //     do {
  //       redo = TS->TimeStep(h, hnext, u);
  //     } while (redo);

  //     u->Print(std::cout);

  //     TS->CommitSolution(h, u);

  //     h = hnext;
  //     i++;

  //     tlast=TS->time();
  //   } while (tout > tlast);

  //   // compute the error with the exact solution
  //   u_ex->putScalar(-1.0/3.0);
  //   u->update(1.0, *u_ex, -1.0);

  //   double norm;
  //   norm = u->normInf();

  //   TS->ReportStatistics_();
  //   CHECK_CLOSE(0.0,norm,1e-3);
  // }
}
