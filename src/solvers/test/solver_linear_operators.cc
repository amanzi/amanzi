#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"

#include "exceptions.hh"
#include "PCG_Operator.hh"

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
    std::cout << "PCG SOLVER..." << std::endl;

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
    std::cout << "GMRES SOLVER..." << std::endl;

    Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
    Epetra_Map* map = new Epetra_Map(2, 0, *comm);

    // create the pcg operator
    Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix());
    AmanziSolvers::PCG_Operator<Matrix, Epetra_Vector, Epetra_Map> gmres(m);

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
}




