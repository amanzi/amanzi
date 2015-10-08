#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"

#include "exceptions.hh"
#include "LinearOperatorFactory.hh"
#include "LinearOperatorPCG.hh"
#include "LinearOperatorGMRES.hh"
#include "LinearOperatorNKA.hh"

using namespace Amanzi;

SUITE(SOLVERS) {
class Matrix {
 public:
  Matrix() {}
  Matrix(const Teuchos::RCP<Epetra_Map>& map) : map_(map) {};
  ~Matrix() {};
  Matrix(const Matrix& other) :
      map_(other.map_) {}
    
  Teuchos::RCP<Matrix> Clone() const {
    return Teuchos::rcp(new Matrix(*this));
  }
    
  virtual int Apply(const Epetra_Vector& v, Epetra_Vector& mv) const { 
    for (int i = 0; i < 5; i++) mv[i] = 2 * v[i];
    for (int i = 1; i < 5; i++) mv[i] -= v[i - 1];
    for (int i = 0; i < 4; i++) mv[i] -= v[i + 1];
    return 0;
  }
  virtual int ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const {
    for (int i = 0; i < 5; i++) hv[i] = v[i];
    return 0;
  }

  virtual const Epetra_Map& DomainMap() const {
    return *map_;
  }
  virtual const Epetra_Map& RangeMap() const {
    return *map_;
  }

 private:
  Teuchos::RCP<Epetra_Map> map_;
};

TEST(PCG_SOLVER) {
  std::cout << "Checking PCG solver..." << std::endl;

  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(5, 0, *comm));

  // create the pcg operator
  Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
  Teuchos::RCP<Matrix> h = m;
  AmanziSolvers::LinearOperatorPCG<Matrix, Epetra_Vector, Epetra_Map> pcg(m, h);
  pcg.Init();

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

  delete comm;
};

TEST(GMRES_SOLVER_LEFT_PRECONDITIONER) {
  std::cout << "\nChecking GMRES solver with LEFT preconditioner..." << std::endl;

  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(5, 0, *comm));

  // create the gmres operator
  for (int loop = 0; loop < 2; loop++) {
    Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
    AmanziSolvers::LinearOperatorGMRES<Matrix, Epetra_Vector, Epetra_Map> gmres(m, m);
    gmres.Init();
    gmres.set_krylov_dim(3 + loop);

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
  }

  delete comm;
};

TEST(GMRES_SOLVER_RIGHT_PRECONDITIONER) {
  std::cout << "\nChecking GMRES solver with RIGHT preconditioner..." << std::endl;

  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(5, 0, *comm));

  Teuchos::ParameterList plist;
  plist.set<std::string>("preconditioning strategy", "right");

  // create the gmres operator
  for (int loop = 0; loop < 2; loop++) {
    Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
    AmanziSolvers::LinearOperatorGMRES<Matrix, Epetra_Vector, Epetra_Map> gmres(m, m);
    gmres.Init(plist);
    gmres.set_krylov_dim(3 + loop);

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
  }

  delete comm;
};

TEST(NKA_SOLVER) {
  std::cout << "Checking NKA solver..." << std::endl;

  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(5, 0, *comm));

  // create the pcg operator
  Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
  AmanziSolvers::LinearOperatorNKA<Matrix, Epetra_Vector, Epetra_Map> nka(m, m);
  nka.Init();

  // initial guess
  Epetra_Vector u(*map);
  u[0] = -1.0;
  u[1] =  1.0;

  // solve
  Epetra_Vector v(*map);
  nka.ApplyInverse(u, v);

  CHECK_CLOSE(-0.1666666666e+0, v[0], 1e-6);
  CHECK_CLOSE( 0.6666666666e+0, v[1], 1e-6);
  CHECK_CLOSE( 0.5e+0, v[2], 1e-6);
  CHECK_CLOSE( 0.3333333333e+0, v[3], 1e-6);
  CHECK_CLOSE( 0.1666666666e+0, v[4], 1e-6);

  delete comm;
};

TEST(SOLVER_FACTORY) {
  std::cout << "Checking solver factory..." << std::endl;

  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(5, 0, *comm));

  Teuchos::ParameterList plist;
  Teuchos::ParameterList& slist = plist.sublist("pcg");
  slist.set<std::string>("iterative method", "pcg");

  // create the pcg operator
  Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
  AmanziSolvers::LinearOperatorFactory<Matrix, Epetra_Vector, Epetra_Map> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Matrix, Epetra_Vector, Epetra_Map> > 
      solver = factory.Create("pcg", plist, m);
  solver->Init();

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

  delete comm;
};

TEST(VERBOSITY_OBJECT) {
  std::cout << "Checking verbosity object..." << std::endl;

  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(5, 0, *comm));

  Teuchos::ParameterList plist;
  Teuchos::ParameterList& slist = plist.sublist("gmres");
  slist.set<std::string>("iterative method", "gmres");
  slist.sublist("gmres parameters").set("size of Krylov space", 5);
  Teuchos::ParameterList& vlist = slist.sublist("VerboseObject");
  vlist.set("Verbosity Level", "extreme");

  // create the pcg operator
  Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
  AmanziSolvers::LinearOperatorFactory<Matrix, Epetra_Vector, Epetra_Map> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Matrix, Epetra_Vector, Epetra_Map> > 
      solver = factory.Create("gmres", plist, m);
  solver->Init();

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

  delete comm;
};

}  // suite SOLVERS




