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

#include "Solver.hh"
#include "SolverFactory.hh"

#include "dbc.hh"
#include "exceptions.hh"
#include "errors.hh"
#include "FnBaseDefs.hh"

#include "ode_fnbase.hh"

using namespace Amanzi;
using namespace Amanzi::AmanziSolvers;


SUITE(MR_GARK_Tests) {
  // data structures for testing
  struct test_data {
    Epetra_MpiComm *comm;

    Teuchos::RCP<Epetra_Vector> init;
    Teuchos::RCP<Epetra_Vector> u;
    Teuchos::RCP<Epetra_Vector> u_dot;
    Teuchos::RCP<Epetra_Vector> u_eval;
    Teuchos::RCP<Epetra_CrsMatrix> A_fast;
    Teuchos::RCP<Epetra_CrsMatrix> A_slow;

    Teuchos::ParameterList plist;

    test_data() {
      std::cout << "Intializing Data" << std::endl;
      // comm and mesh for maps
      comm = new Epetra_MpiComm(MPI_COMM_SELF);
      Epetra_Map map(2,0,*comm);
      Epetra_Map map_matrix(2, 0, *comm);
      init = Teuchos::rcp(new Epetra_Vector(map));

      // u, u_dot, and exact soln
      u = Teuchos::rcp(new Epetra_Vector(*init));
      u_dot = Teuchos::rcp(new Epetra_Vector(*init));
      u_eval = Teuchos::rcp(new Epetra_Vector(*init));

      A_fast = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map_matrix, 2, true));
      A_slow = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map_matrix, 2, true));

      
      std::cout << "Intialized Data" << std::endl;
      
      const int row_index[2] = {0,1};
      double values[2] = {2.0,1.0};
      A_fast->InsertGlobalValues(0, 2, values, row_index);

      values[0] = 0.0;
      values[1] = 2.0;
      A_fast->InsertGlobalValues(1, 2, values, row_index);

      
      values[0] = 1.0;
      values[1] = 0.0;
      A_slow->InsertGlobalValues(0, 2, values, row_index);

      
      values[0] = 1.0;
      values[1] = 1.0;
      A_slow->InsertGlobalValues(1, 2, values, row_index);

      A_fast->FillComplete();
      A_slow->FillComplete();

      std::cout << "Set Matrices Data" << std::endl;

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
    
    MRG_EXIM_linear2d_ODE ToyProblem(1., 1., true, A_fast, A_slow, comm);

    // create the time stepper
    Teuchos::RCP<Amanzi::MRG_EXIM_TI<Epetra_MultiVector, Epetra_BlockMap> > TS =
        Teuchos::rcp(new MRG_EXIM_TI<Epetra_MultiVector, Epetra_BlockMap>(ToyProblem, MrGark_EXIM_2, plist, init ));

    // initial value
    u->PutScalar(1.0);
    u_dot->PutScalar(1.0);

    // initial time
    double t=0.0;

    // final time
    double tout = 2.0;

    // initial time step
    double h = 1.0e-3;
    double hnext;

    //TimeStep Ratio
    int M  = 3;

    // iterate until the final time
    int i=0;
    double tlast = t;

    std::cout << "starting time integration" << std::endl;
    do {
      if (tlast + h > tout) {
        std::cout << "adjusting h, to hit the final time exactly...\n";
        h = tout - tlast;
      }

      int redo = 0;
      do {
        redo = TS->TimeStep(tlast, h, M, u, u_eval);
      } while (redo);

      u_eval->Print(std::cout);

      h = hnext;
      i++;

      *u = *u_eval;

    } while (tout > tlast);

    ToyProblem.exact_rhs(tout, u_eval);

    u_eval->Print(std::cout);
    u->Print(std::cout);

    u_eval->Update(-1.0, *u, 1.0);

    double norm;
    u_eval->NormInf(&norm);

    std::cout << "Absolute Error: " << norm << std::endl;

    CHECK_CLOSE(0.0,norm,1e-3);
  }

 
}


