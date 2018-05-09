#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"

#include "exceptions.hh"
#include "LinearOperatorFactory.hh"
#include "LinearOperatorPCG.hh"
#include "LinearOperatorGMRES.hh"
#include "LinearOperatorNKA.hh"
#include "LinearOperatorBelosGMRES.hh"
#include "LinearOperatorAmesos.hh"

using namespace Amanzi;

SUITE(SOLVERS) {
class Matrix {
 public:
  Matrix() {};
  Matrix(const Teuchos::RCP<Epetra_Map>& map) : map_(map) {
    x_[0] = 0.00699270335645641;
    x_[1] = 0.01398540671291281;
    x_[2] = 0.02079636439636044;
    x_[3] = 0.02688033938777295;
    x_[4] = 0.03122225970045909;
  };
  ~Matrix() {};
  Matrix(const Matrix& other) : map_(other.map_) {};
    
  Teuchos::RCP<Matrix> Clone() const {
    return Teuchos::rcp(new Matrix(*this));
  }
    
  // 5-point FD stensil
  virtual int Apply(const Epetra_MultiVector& v, Epetra_MultiVector& mv) const { 
    int n = std::pow(v.Map().NumMyElements(), 0.5);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        int k = j * n + i;
        mv[0][k] = 4 * v[0][k];
        if (i > 0) mv[0][k] -= v[0][k - 1];
        if (j > 0) mv[0][k] -= v[0][k - n];

        if (i < n - 1) mv[0][k] -= v[0][k + 1];
        if (j < n - 1) mv[0][k] -= v[0][k + n];
      }
    }
    return 0;
  }

  virtual int ApplyInverse(const Epetra_MultiVector& v, Epetra_MultiVector& hv) const {
    int n = v.Map().NumMyElements();
    for (int i = 0; i < n; i++) hv[0][i] = v[0][i];
    return 0;
  }

  // 3-point FD stensil
  void Init() {
    int n = map_->NumMyElements();
    A_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *map_, *map_, 3));
    for (int i = 0; i < n; i++) {
      int indices[3];
      double values[3] = {double(-i), double(2 * i + 1), double(-i - 1)};
      for (int k = 0; k < 3; k++) indices[k] = i + k - 1; 
      A_->InsertMyValues(i, 3, values, indices);
    }
    A_->FillComplete(*map_, *map_);
  }

  Teuchos::RCP<Epetra_CrsMatrix> A() const { return A_; }

  virtual const Epetra_Map& DomainMap() const { return *map_; }
  virtual const Epetra_Map& RangeMap() const { return *map_; }
  virtual double* x() { return x_; }

 private:
  Teuchos::RCP<Epetra_Map> map_;
  double x_[5];
  Teuchos::RCP<Epetra_CrsMatrix> A_;
};


TEST(PCG_SOLVER) {
  std::cout << "Checking PCG solver..." << std::endl;

  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(100, 0, *comm));

  // create the pcg operator
  Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
  Teuchos::RCP<Matrix> h = m;
  AmanziSolvers::LinearOperatorPCG<Matrix, Epetra_Vector, Epetra_Map> pcg(m, h);
  pcg.Init();

  // initial guess
  Epetra_Vector u(*map);
  u[55] = 1.0;

  // solve
  Epetra_Vector v(*map);
  int ierr = pcg.ApplyInverse(u, v);
  CHECK(ierr > 0);

  for (int i = 0; i < 5; i++) CHECK_CLOSE((m->x())[i], v[i], 1e-6);

  delete comm;
};

TEST(GMRES_SOLVER_LEFT_PRECONDITIONER) {
  std::cout << "\nChecking GMRES solver with LEFT preconditioner..." << std::endl;

  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(100, 0, *comm));

  // create the gmres operator
  for (int loop = 0; loop < 2; loop++) {
    Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
    AmanziSolvers::LinearOperatorGMRES<Matrix, Epetra_Vector, Epetra_Map> gmres(m, m);
    gmres.Init();
    gmres.set_krylov_dim(15 + loop * 5);
    gmres.set_tolerance(1e-12);

    // initial guess
    Epetra_Vector u(*map);
    u[55] = 1.0;

    // solve
    Epetra_Vector v(*map);
    int ierr = gmres.ApplyInverse(u, v);
    CHECK(ierr > 0);
    for (int i = 0; i < 5; i++) CHECK_CLOSE((m->x())[i], v[i], 1e-6);
  }

  delete comm;
};

TEST(GMRES_SOLVER_RIGHT_PRECONDITIONER) {
  std::cout << "\nChecking GMRES solver with RIGHT preconditioner..." << std::endl;

  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(100, 0, *comm));

  Teuchos::ParameterList plist;
  plist.set<std::string>("preconditioning strategy", "right");

  // create the gmres operator
  for (int loop = 0; loop < 2; loop++) {
    Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
    AmanziSolvers::LinearOperatorGMRES<Matrix, Epetra_Vector, Epetra_Map> gmres(m, m);
    gmres.Init(plist);
    gmres.set_krylov_dim(15 + loop * 5);
    gmres.set_tolerance(1e-12);

    // initial guess
    Epetra_Vector u(*map);
    u[55] = 1.0;

    // solve
    Epetra_Vector v(*map);
    int ierr = gmres.ApplyInverse(u, v);
    CHECK(ierr > 0);

    for (int i = 0; i < 5; i++) CHECK_CLOSE((m->x())[i], v[i], 1e-6);
  }

  delete comm;
};

TEST(GMRES_SOLVER_DEFLATION) {
  std::cout << "\nChecking GMRES solver with deflated restart..." << std::endl;

  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(100, 0, *comm));

  Teuchos::ParameterList plist;
  plist.set<int>("maximum size of deflation space", 5);
  plist.set<int>("maximum number of iterations", 200);
  plist.set<int>("size of Krylov space", 15);
  plist.set<double>("error tolerance", 1e-12);
  Teuchos::ParameterList& vlist = plist.sublist("verbose object");
  vlist.set("verbosity level", "extreme");

  // create the gmres operator
  Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
  AmanziSolvers::LinearOperatorGMRES<Matrix, Epetra_Vector, Epetra_Map> gmres(m, m);
  gmres.Init(plist);

  // initial guess
  Epetra_Vector u(*map);
  u[55] = 1.0;

  // solve
  Epetra_Vector v(*map);
  int ierr = gmres.ApplyInverse(u, v);
  CHECK(ierr > 0);

  for (int i = 0; i < 5; i++) CHECK_CLOSE((m->x())[i], v[i], 1e-6);

  delete comm;
};

TEST(NKA_SOLVER) {
  std::cout << "\nChecking NKA solver..." << std::endl;

  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(100, 0, *comm));

  // create the pcg operator
  Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
  AmanziSolvers::LinearOperatorNKA<Matrix, Epetra_Vector, Epetra_Map> nka(m, m);
  nka.Init();

  // initial guess
  Epetra_Vector u(*map);
  u[55] = 1.0;

  // solve
  Epetra_Vector v(*map);
  int ierr = nka.ApplyInverse(u, v);
  CHECK(ierr > 0);

  for (int i = 0; i < 5; i++) CHECK_CLOSE((m->x())[i], v[i], 1e-6);

  delete comm;
};

TEST(BELOS_GMRES_SOLVER) {
  std::cout << "\nChecking Belos GMRES solver..." << std::endl;

  Teuchos::ParameterList plist;
  Teuchos::ParameterList& vo = plist.sublist("VerboseObject");
  vo.set("Verbosity Level", "high");
  plist.set<int>("size of Krylov space", 15);
  plist.set<double>("error tolerance", 1e-12);
  
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(100, 0, *comm));

  // create the operator
  Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
  AmanziSolvers::LinearOperatorBelosGMRES<Matrix, Epetra_Vector, Epetra_Map> gmres(m, m);
  gmres.Init(plist);

  // initial guess
  Epetra_Vector u(*map);
  u[55] = 1.0;

  // solve
  Epetra_Vector v(*map);
  int ierr = gmres.ApplyInverse(u, v);
  CHECK(ierr > 0);

  for (int i = 0; i < 5; i++) CHECK_CLOSE((m->x())[i], v[i], 1e-6);

  delete comm;
};

TEST(AMESOS_SOLVER) {
  std::cout << "\nChecking Amesos solver..." << std::endl;

  Teuchos::ParameterList plist;
  Teuchos::ParameterList& vo = plist.sublist("VerboseObject");
  vo.set("Verbosity Level", "high");
  plist.set<std::string>("solver name", "Amesos_Klu");
  
  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(10, 0, *comm));

  // create the operator
  Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
  m->Init();
  AmanziSolvers::LinearOperatorAmesos<Matrix, Epetra_MultiVector, Epetra_Map> superlu(m, m);
  superlu.Init(plist);

  // initial guess
  Epetra_Vector u(*map);
  u[5] = 1.0;

  // solve 1D problem
  Epetra_Vector v(*map);
  int ierr = superlu.ApplyInverse(u, v);
  CHECK(ierr > 0);

  double residual = 11 * v[5] - 5 * v[4] - 6 * v[6];
  CHECK_CLOSE(residual, 1.0, 1e-12); 

  delete comm;
};

TEST(SOLVER_FACTORY) {
  std::cout << "\nChecking solver factory..." << std::endl;

  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(100, 0, *comm));

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
  u[55] = 1.0;

  // solve
  Epetra_Vector v(*map);
  int ierr = solver->ApplyInverse(u, v);
  CHECK(ierr > 0);

  for (int i = 0; i < 5; i++) CHECK_CLOSE((m->x())[i], v[i], 1e-6);

  delete comm;
};

TEST(VERBOSITY_OBJECT) {
  std::cout << "\nChecking verbosity object..." << std::endl;

  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(100, 0, *comm));

  Teuchos::ParameterList plist;
  Teuchos::ParameterList& slist = plist.sublist("gmres");
  slist.set<std::string>("iterative method", "gmres");
  slist.sublist("gmres parameters")
      .set("size of Krylov space", 50)
      .set<double>("error tolerance", 1e-12);
  slist.sublist("gmres parameters")
       .sublist("verbose object").set("verbosity level", "extreme");

  // create the pcg operator
  Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
  AmanziSolvers::LinearOperatorFactory<Matrix, Epetra_Vector, Epetra_Map> factory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Matrix, Epetra_Vector, Epetra_Map> > 
      solver = factory.Create("gmres", plist, m);

  // initial guess
  Epetra_Vector u(*map);
  u[55] = 1.0;

  // solve
  Epetra_Vector v(*map);
  int ierr = solver->ApplyInverse(u, v);
  CHECK(ierr > 0);

  for (int i = 0; i < 5; i++) CHECK_CLOSE((m->x())[i], v[i], 1e-6);

  delete comm;
};

}  // suite SOLVERS




