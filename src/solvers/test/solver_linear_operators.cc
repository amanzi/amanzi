#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"

#include "exceptions.hh"
#include "LinearOperatorFactory.hh"
#include "LinearOperatorPCG.hh"
#include "LinearOperatorGMRES.hh"

using namespace Amanzi;

SUITE(SOLVERS) {
  class Matrix {
   public:
    Matrix() {};
    ~Matrix() {};

    void Apply(const Epetra_Vector& v, Epetra_Vector& mv) const { 
      for (int i = 0; i < 5; i++) mv[i] = 2 * v[i];
      for (int i = 1; i < 5; i++) mv[i] -= v[i - 1];
      for (int i = 0; i < 4; i++) mv[i] -= v[i + 1];
    }
    void ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const {
      for (int i = 0; i < 5; i++) hv[i] = v[i];
    }
  };

  TEST(PCG_SOLVER) {
    std::cout << "Checking PCG solver..." << std::endl;

    Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
    Epetra_Map* map = new Epetra_Map(5, 0, *comm);

    // create the pcg operator
    Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix());
    AmanziSolvers::LinearOperatorPCG<Matrix, Epetra_Vector, Epetra_Map> pcg(m, m);

    // initial guess
    Epetra_Vector u(*map);
    u[0] = -1.0;
    u[1] =  1.0;

    // solve
    Epetra_Vector v(*map);
    pcg.ApplyInverse(u, v);

    CHECK_CLOSE(-0.1666666666e+0, v[0], 1e-6);
    CHECK_CLOSE( 0.6666666666e+0, v[1], 1e-6);
    CHECK_CLOSE( 0.5e+0, v[2], 1e-6);
    CHECK_CLOSE( 0.3333333333e+0, v[3], 1e-6);
    CHECK_CLOSE( 0.1666666666e+0, v[4], 1e-6);
  };

  TEST(GMRES_SOLVER) {
    std::cout << "Checking GMRES solver..." << std::endl;

    Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
    Epetra_Map* map = new Epetra_Map(5, 0, *comm);

    // create the pcg operator
    Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix());
    AmanziSolvers::LinearOperatorGMRES<Matrix, Epetra_Vector, Epetra_Map> gmres(m, m);

    // initial guess
    Epetra_Vector u(*map);
    u[0] = -1.0;
    u[1] =  1.0;

    // solve
    Epetra_Vector v(*map);
    gmres.ApplyInverse(u, v);

    CHECK_CLOSE(-0.1666666666e+0, v[0], 1e-6);
    CHECK_CLOSE( 0.6666666666e+0, v[1], 1e-6);
    CHECK_CLOSE( 0.5e+0, v[2], 1e-6);
    CHECK_CLOSE( 0.3333333333e+0, v[3], 1e-6);
    CHECK_CLOSE( 0.1666666666e+0, v[4], 1e-6);
  };

  TEST(SOLVER_FACTORY) {
    std::cout << "Checking solver factory..." << std::endl;

    Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
    Epetra_Map* map = new Epetra_Map(5, 0, *comm);

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

    CHECK_CLOSE(-0.1666666666e+0, v[0], 1e-6);
    CHECK_CLOSE( 0.6666666666e+0, v[1], 1e-6);
    CHECK_CLOSE( 0.5e+0, v[2], 1e-6);
    CHECK_CLOSE( 0.3333333333e+0, v[3], 1e-6);
    CHECK_CLOSE( 0.1666666666e+0, v[4], 1e-6);
  };

  TEST(VERBOSITY_OBJECT) {
    std::cout << "Checking verbosity object..." << std::endl;

    Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
    Epetra_Map* map = new Epetra_Map(5, 0, *comm);

    Teuchos::ParameterList plist;
    Teuchos::ParameterList& slist = plist.sublist("gmres");
    slist.set<string>("iterative method", "gmres");
    Teuchos::ParameterList& vlist = slist.sublist("VerboseObject");
    vlist.set("Verbosity Level", "extreme");

    // create the pcg operator
    Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix());
    AmanziSolvers::LinearOperatorFactory<Matrix, Epetra_Vector, Epetra_Map> factory;
    Teuchos::RCP<AmanziSolvers::LinearOperator<Matrix, Epetra_Vector, Epetra_Map> > 
        solver = factory.Create("gmres", plist, m);

    // initial guess
    Epetra_Vector u(*map);
    u[0] = -1.0;
    u[1] =  1.0;

    // solve
    Epetra_Vector v(*map);
    solver->ApplyInverse(u, v);

    CHECK_CLOSE(-0.1666666666e+0, v[0], 1e-6);
    CHECK_CLOSE( 0.6666666666e+0, v[1], 1e-6);
    CHECK_CLOSE( 0.5e+0, v[2], 1e-6);
    CHECK_CLOSE( 0.3333333333e+0, v[3], 1e-6);
    CHECK_CLOSE( 0.1666666666e+0, v[4], 1e-6);
  };

}




