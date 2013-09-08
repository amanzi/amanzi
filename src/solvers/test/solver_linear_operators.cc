#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"

#include "exceptions.hh"
#include "LinearOperatorFactory.hh"
#include "PCG_Operator.hh"
#include "GMRES_Operator.hh"

using namespace Amanzi;

SUITE(SOLVERS) {
  class Matrix {
   public:
    Matrix() {};
    ~Matrix() {};

    void Apply(const Epetra_Vector& v, Epetra_Vector& mv) const { 
      mv[0] = 2 * v[0] - v[1]; 
      mv[1] = 2 * v[1] - v[0]; 
    }
    void ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const {
      hv[0] = v[0]; 
      hv[1] = v[1]; 
    }
  };

  TEST(PCG_SOLVER) {
    std::cout << "Checking PCG solver..." << std::endl;

    Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
    Epetra_Map* map = new Epetra_Map(2, 0, *comm);

    // create the pcg operator
    Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix());
    AmanziSolvers::PCG_Operator<Matrix, Epetra_Vector, Epetra_Map> pcg(m);

    // initial guess
    Epetra_Vector u(*map);
    u[0] = -1.0;
    u[1] =  1.0;

    // solve
    Epetra_Vector v(*map);
    pcg.ApplyInverse(u, v);
    CHECK_CLOSE(-0.33333333e+0, v[0], 1e-6);
    CHECK_CLOSE( 0.33333333e+0, v[1], 1e-6);
  };

  TEST(GMRES_SOLVER) {
    std::cout << "Checking GMRES solver..." << std::endl;

    Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
    Epetra_Map* map = new Epetra_Map(2, 0, *comm);

    // create the pcg operator
    Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix());
    AmanziSolvers::GMRES_Operator<Matrix, Epetra_Vector, Epetra_Map> gmres(m);

    // initial guess
    Epetra_Vector u(*map);
    u[0] = -1.0;
    u[1] =  1.0;

    // solve
    Epetra_Vector v(*map);
    gmres.ApplyInverse(u, v);
    CHECK_CLOSE(-0.33333333e+0, v[0], 1e-6);
    CHECK_CLOSE( 0.33333333e+0, v[1], 1e-6);
  };

  TEST(SOLVER_FACTORY) {
    std::cout << "Checking solver factory..." << std::endl;

    Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
    Epetra_Map* map = new Epetra_Map(2, 0, *comm);

    Teuchos::ParameterList plist;
    Teuchos::ParameterList& slist = plist.sublist("pcg");
    slist.set<string>("iterative method", "pcg");

    // create the pcg operator
    Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix());
    AmanziSolvers::LinearOperatorFactory<Matrix, Epetra_Vector, Epetra_Map> factory;
    Teuchos::RCP<AmanziSolvers::LinearOperator<Matrix, Epetra_Vector, Epetra_Map> > 
        solver = factory.Create("pcg", plist, m);

    // initial guess
    Epetra_Vector u(*map);
    u[0] = -1.0;
    u[1] =  1.0;

    // solve
    Epetra_Vector v(*map);
    solver->ApplyInverse(u, v);
    CHECK_CLOSE(-0.33333333e+0, v[0], 1e-6);
    CHECK_CLOSE( 0.33333333e+0, v[1], 1e-6);
  };

  TEST(VERBOSITY_OBJECT) {
    std::cout << "Checking verbosity object..." << std::endl;

    Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
    Epetra_Map* map = new Epetra_Map(2, 0, *comm);

    Teuchos::ParameterList plist;
    Teuchos::ParameterList& slist = plist.sublist("pcg");
    slist.set<string>("iterative method", "pcg");
    Teuchos::ParameterList& vlist = slist.sublist("VerboseObject");
    vlist.set("Verbosity Level", "extreme");

    // create the pcg operator
    Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix());
    AmanziSolvers::LinearOperatorFactory<Matrix, Epetra_Vector, Epetra_Map> factory;
    Teuchos::RCP<AmanziSolvers::LinearOperator<Matrix, Epetra_Vector, Epetra_Map> > 
        solver = factory.Create("pcg", plist, m);

    // initial guess
    Epetra_Vector u(*map);
    u[0] = -1.0;
    u[1] =  1.0;

    // solve
    Epetra_Vector v(*map);
    solver->ApplyInverse(u, v);
    CHECK_CLOSE(-0.33333333e+0, v[0], 1e-6);
    CHECK_CLOSE( 0.33333333e+0, v[1], 1e-6);
  };

}




