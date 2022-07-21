#include <ios>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerbosityLevel.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"

#include "MRG_EXIM_FnBase.hh"
#include "MRG_EXIM_SolverFnBase.hh"
#include "MRG_EXIM_TI.hh"
#include "MRG_EXIM_ode_fnbase.hh"
#include "PartitionFnBase.hh"

#include "MRG_IMIM_FnBase.hh"
#include "MRG_IMIM_ode_fnbase.hh"
#include "MRG_IMIM_SolverFnBase_Fast.hh"
#include "MRG_IMIM_SolverFnBase_Full.hh"
#include "MRG_IMIM_TI.hh"

#include "dbc.hh"
#include "exceptions.hh"
#include "errors.hh"
#include "FnBaseDefs.hh"


using namespace Amanzi;
using namespace AmanziSolvers;


SUITE(MR_GARK_Tests) {
  // data structures for testing
  struct test_data {
    Epetra_MpiComm *comm;

    Teuchos::RCP<Epetra_Vector> init;
    Teuchos::RCP<Epetra_MultiVector> init_full;
    Teuchos::RCP<Epetra_Vector> u;
    Teuchos::RCP<Epetra_Vector> u_eval;
    Teuchos::RCP<Epetra_CrsMatrix> A_fast;
    Teuchos::RCP<Epetra_CrsMatrix> A_slow;

    Teuchos::ParameterList plist;

    test_data() {
      // comm and mesh for maps
      comm = new Epetra_MpiComm(MPI_COMM_SELF);
      Epetra_Map map(2,0,*comm);
      Epetra_Map map_matrix(2, 0, *comm);
      Epetra_Map map_full(4, 0, *comm);
      init = Teuchos::rcp(new Epetra_Vector(map));
      init_full = Teuchos::rcp(new Epetra_MultiVector(map_full, 2));

      // u, u_dot, and exact soln
      u = Teuchos::rcp(new Epetra_Vector(*init));
      u_eval = Teuchos::rcp(new Epetra_Vector(*init));

      A_fast = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map_matrix, 2, true));
      A_slow = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map_matrix, 2, true));
      
      const int row_index[2] = {0,1};
      double values[2] = {-3.0,1.0};
      A_fast->InsertGlobalValues(0, 2, values, row_index);

      values[0] = 0.0;
      values[1] = 1.0;
      A_fast->InsertGlobalValues(1, 2, values, row_index);

      
      values[0] = 1.0;
      values[1] = 0.0;
      A_slow->InsertGlobalValues(0, 2, values, row_index);

      
      values[0] = 2.0;
      values[1] = -10.0;
      A_slow->InsertGlobalValues(1, 2, values, row_index);

      A_fast->FillComplete();
      A_slow->FillComplete();

    }
    ~test_data() {
      delete comm;
    }
  };


  TEST_FIXTURE(test_data, ODE2D_MRG_EXIM_NKA_TrueJacobian) {
    std::cout << "Test: ODE 2D on MR-GARK EXIM (Order 2)" << std::endl;
    std::cout << "\twith NKA, PC is True Jacobian (Fixed Step)" << std::endl;

    // set the parameter list for BDF1
    // set the parameter list for BDF1
    // set up solver params
    plist.set("solver type", "nka");
    plist.sublist("nka parameters").set("limit iterations", 20);
    plist.sublist("nka parameters").set("nonlinear tolerance",1e-10);
    plist.sublist("nka parameters").set("diverged tolerance",1.0e4);
    plist.sublist("nka parameters").set("convergence monitor","monitor update");
    // create the PDE problem

    MRG_EXIM_linear2d_ODE ToyProblem(0.0, 1.0, true, A_fast, A_slow, comm);

    std::cout << "Slow Partition" << std::endl;
    A_slow->Print(std::cout);
    std::cout << "Fast Partition" << std::endl;
    A_fast->Print(std::cout);

    // create the time stepper
    Teuchos::RCP<Amanzi::MRG_EXIM_TI<Epetra_Vector, Epetra_BlockMap> > TS =
        Teuchos::rcp(new MRG_EXIM_TI<Epetra_Vector, Epetra_BlockMap>(ToyProblem, MrGark_EXIM_2, plist, init ));

    // initial value
    u->PutScalar(1.0);

    // initial time
    double t=0.0;

    // final time
    double tout = 2.0;

    // initial time step
    double h = 1.0e-3;
    double hnext;

    //TimeStep Ratio
    int M  = 4;

    // iterate until the final time
    int i=0;
    double tlast = t;

    std::cout << "starting time integration" << std::endl;
    std::cout << "Constant Step Size of : " << h << std::endl;
    std::cout << "M (Stepsize Ratio) : " << M << std::endl;
    do {
      if (tlast + h > tout) {
        std::cout << "adjusting h, to hit the final time exactly...\n";
        h = tout - tlast;
      }

      int redo = 0;
      do {
        redo = TS->TimeStep(tlast, h, M, u, u_eval);
      } while (redo);

      i++;

      *u = *u_eval;

      // u->Print(std::cout);

      tlast += h;
    } while (tout > tlast);

    ToyProblem.exact_rhs(tout, u_eval);

    std::cout << "Exact Solution" << std::endl;
    u_eval->Print(std::cout);

    std::cout << "Numerical Solution" << std::endl;
    u->Print(std::cout);

    u_eval->Update(-1.0, *u, 1.0);

    double norm;
    u_eval->NormInf(&norm);

    std::cout << "Absolute Error: " << norm << std::endl;

    CHECK_CLOSE(0.0,norm,1e-3);
  }

  TEST_FIXTURE(test_data, ODE2D_MRG_IMIM_NKA_TrueJacobian) {
    std::cout << "Test: ODE 2D on MR-GARK IMIM (Order 2)" << std::endl;
    std::cout << "\twith NKA, PC is True Jacobian (Fixed Step)" << std::endl;

    // set the parameter list for BDF1
    // set up solver params
    plist.set("solver type", "nka");
    plist.sublist("nka parameters").set("limit iterations", 20);
    plist.sublist("nka parameters").set("nonlinear tolerance",1e-10);
    plist.sublist("nka parameters").set("diverged tolerance",1.0e4);
    plist.sublist("nka parameters").set("convergence monitor","monitor update");
    // create the PDE problem

    MRG_IMIM_linear2d_ODE ToyProblem(0.0, 1.0, true, A_fast, A_slow, comm);

    std::cout << "Problem Initalized" << std::endl;

    // create the time stepper
    Teuchos::RCP<Amanzi::MRG_IMIM_TI<Epetra_Vector, Epetra_BlockMap, Epetra_MultiVector, Epetra_BlockMap> > TS =
        Teuchos::rcp(new MRG_IMIM_TI<Epetra_Vector, Epetra_BlockMap, Epetra_MultiVector, Epetra_BlockMap>(ToyProblem, MrGark_IMIM_2, plist, plist, init, init_full));

    std::cout << "Time Integrator Intialized" << std::endl;

    // initial value
    u->PutScalar(1.0);

    // initial time
    double t=0.0;

    // final time
    double tout = 2.0;

    // initial time step
    double h = 1.0e-3;
    double hnext;

    //TimeStep Ratio
    int M  = 1;

    // iterate until the final time
    int i=0;
    double tlast = t;

    std::cout << "starting time integration" << std::endl;
    std::cout << "Constant Step Size of : " << h << std::endl;
    std::cout << "M (Stepsize Ratio) : " << M << std::endl;
    do {
      if (tlast + h > tout) {
        std::cout << "adjusting h, to hit the final time exactly...\n";
        h = tout - tlast;
      }

      int redo = 0;
      do {
        redo = TS->TimeStep(tlast, h, M, u, u_eval);
      } while (redo);

      i++;

      *u = *u_eval;

      // u->Print(std::cout);

      tlast += h;
    } while (tout > tlast);

    ToyProblem.exact_rhs(tout, u_eval);

    std::cout << "Exact Solution" << std::endl;
    u_eval->Print(std::cout);

    std::cout << "Numerical Solution" << std::endl;
    u->Print(std::cout);

    u_eval->Update(-1.0, *u, 1.0);

    double norm;
    u_eval->NormInf(&norm);

    std::cout << "Absolute Error: " << norm << std::endl;

    CHECK_CLOSE(0.0,norm,1e-3);
  }

 
}


