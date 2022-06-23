#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"

#include "exceptions.hh"
#include "SuperMap.hh"
#include "InverseFactory.hh"
#include "IterativeMethodPCG.hh"
#include "IterativeMethodGMRES.hh"
#include "IterativeMethodNKA.hh"
#include "IterativeMethodBelos.hh"
#include "DirectMethodAmesos.hh"
#include "DirectMethodAmesos2.hh"

using namespace Amanzi;

SUITE(SOLVERS) {
class Matrix {
 public:
  using Vector_t = Epetra_Vector;
  using VectorSpace_t = Epetra_Map;

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

  // 5-point FD stencil
  int Apply(const Epetra_Vector& v, Epetra_Vector& mv) const {
    int n = std::pow(v.Map().NumMyElements(), 0.5);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        int k = j * n + i;
        mv[k] = 4 * v[k];
        if (i > 0) mv[k] -= v[k - 1];
        if (j > 0) mv[k] -= v[k - n];

        if (i < n - 1) mv[k] -= v[k + 1];
        if (j < n - 1) mv[k] -= v[k + n];
      }
    }
    return 0;
  }

  void ComputeInverse() {}
  void InitializeInverse() {}

  int ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const {
    int n = v.Map().NumMyElements();
    for (int i = 0; i < n; i++) hv[i] = v[i];
    return 0;
  }

  // 3-point FD stencil (for direct solvers)
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

  int SymbolicAssembleMatrix() { return 0; }
  int AssembleMatrix() { return 0;}

  // partial consistency with Operators'interface
  Teuchos::RCP<Epetra_CrsMatrix> A() { return A_; }
  Teuchos::RCP<Amanzi::Operators::SuperMap> get_supermap() const { return Teuchos::null; }

  const Epetra_Map& DomainMap() const { return *map_; }
  const Epetra_Map& RangeMap() const { return *map_; }
  double* x() { return x_; }

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
  AmanziSolvers::IterativeMethodPCG<Matrix, Matrix, Epetra_Vector, Epetra_Map> pcg;
  pcg.set_matrices(m,h);
  Teuchos::ParameterList plist;
  pcg.set_inverse_parameters(plist);

  // initial guess
  Epetra_Vector u(*map);
  u[55] = 1.0;

  // solve
  Epetra_Vector v(*map);
  int ierr = pcg.ApplyInverse(u, v);
  CHECK(ierr == 0);

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
    AmanziSolvers::IterativeMethodGMRES<Matrix,Matrix,Epetra_Vector,Epetra_Map> gmres;
    Teuchos::ParameterList plist;
    gmres.set_inverse_parameters(plist);
    gmres.set_matrices(m,m);
    gmres.set_krylov_dim(15 + loop * 5);
    gmres.set_tolerance(1e-12);

    // initial guess
    Epetra_Vector u(*map);
    u[55] = 1.0;

    // solve
    Epetra_Vector v(*map);
    int ierr = gmres.ApplyInverse(u, v);
    CHECK(ierr == 0);
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
  plist.set<bool>("release Krylov vectors", true);

  // create the gmres operator
  for (int loop = 0; loop < 2; loop++) {
    Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
    AmanziSolvers::IterativeMethodGMRES<Matrix,Matrix,Epetra_Vector,Epetra_Map> gmres;
    gmres.set_inverse_parameters(plist);
    gmres.set_matrices(m,m);
    gmres.set_krylov_dim(15 + loop * 5);
    gmres.set_tolerance(1e-12);

    // initial guess
    Epetra_Vector u(*map);
    u[55] = 1.0;

    // solve
    Epetra_Vector v(*map);
    int ierr = gmres.ApplyInverse(u, v);
    CHECK(ierr == 0);

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
  AmanziSolvers::IterativeMethodGMRES<Matrix,Matrix,Epetra_Vector,Epetra_Map> gmres;
  gmres.set_matrices(m, m);
  gmres.set_inverse_parameters(plist);

  // initial guess
  Epetra_Vector u(*map);
  u[55] = 1.0;

  // solve
  Epetra_Vector v(*map);
  int ierr = gmres.ApplyInverse(u, v);
  CHECK(ierr == 0);

  for (int i = 0; i < 5; i++) CHECK_CLOSE((m->x())[i], v[i], 1e-6);

  delete comm;
};

TEST(NKA_SOLVER) {
  std::cout << "\nChecking NKA solver..." << std::endl;

  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(100, 0, *comm));

  // create NKA operator
  Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
  AmanziSolvers::IterativeMethodNKA<Matrix,Matrix,Epetra_Vector,Epetra_Map> nka;
  nka.set_matrices(m,m);
  Teuchos::ParameterList plist;
  plist.set("error tolerance", 1.e-13);
  plist.set("maximum number of iterations", 200);
  plist.sublist("verbose object").set("verbosity level", "high");
  nka.set_inverse_parameters(plist);

  // initial guess
  Epetra_Vector u(*map);
  u[55] = 1.0;

  // solve
  Epetra_Vector v(*map);
  int ierr = nka.ApplyInverse(u, v);
  CHECK(ierr == 0);

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

  Epetra_MpiComm comm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(100, 0, comm));

  // create the operator
  Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
  AmanziSolvers::IterativeMethodBelos<Matrix,Matrix,Epetra_Vector,Epetra_Map> gmres;
  gmres.set_matrices(m, m);
  gmres.set_inverse_parameters(plist);
  gmres.InitializeInverse();
  gmres.ComputeInverse();

  // initial guess
  Epetra_Vector u(*map);
  u[55] = 1.0;

  // solve
  Epetra_Vector v(*map);
  int ierr = gmres.ApplyInverse(u, v);
  CHECK(ierr == 0);

  for (int i = 0; i < 5; i++) CHECK_CLOSE((m->x())[i], v[i], 1e-6);
};

TEST(AMESOS_SOLVER) {
  std::cout << "\nChecking Amesos solver..." << std::endl;

  Teuchos::ParameterList plist;
  Teuchos::ParameterList& vo = plist.sublist("VerboseObject");
  vo.set("Verbosity Level", "high");
  plist.set<std::string>("solver name", "Amesos_Klu")
       .set<int>("amesos version", 1);

  Epetra_MpiComm comm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(10, 0, comm));

  // create the operator
  Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
  m->Init();

  // initial guess
  Epetra_Vector v(*map), u(*map);
  u[5] = 1.0;

  // Amesos1
  {
    AmanziSolvers::DirectMethodAmesos klu;
    klu.set_inverse_parameters(plist);
    klu.set_matrix(m->A());
    klu.InitializeInverse();
    klu.ComputeInverse();

    int ierr = klu.ApplyInverse(u, v);
    CHECK(ierr == 0);

    double residual = 11 * v[5] - 5 * v[4] - 6 * v[6];
    CHECK_CLOSE(1.0, residual, 1e-12);
  }

  // Amesos2
  {
    AmanziSolvers::DirectMethodAmesos2 klu;
    plist.set<std::string>("solver name", "Klu2")
      .set<int>("amesos version", 2)
      .set("method", "amesos2: klu");
    klu.set_inverse_parameters(plist);
    klu.set_matrix(m->A());
    klu.InitializeInverse();
    klu.ComputeInverse();

    int ierr = klu.ApplyInverse(u, v);
    CHECK(ierr == 0);

    double residual = 11 * v[5] - 5 * v[4] - 6 * v[6];
    CHECK_CLOSE(1.0, residual, 1e-12);
  }
};

TEST(SOLVER_FACTORY_NO_PC) {
  std::cout << "\nChecking solver factory..." << std::endl;

  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(100, 0, *comm));

  Teuchos::ParameterList plist;
  Teuchos::ParameterList& slist = plist.sublist("pcg");
  slist.set<std::string>("iterative method", "pcg");

  // create the pcg operator
  Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
  auto solver = AmanziSolvers::createInverse<Matrix,Epetra_Vector,Epetra_Map>(
      "pcg", plist, m);

  // initial guess
  Epetra_Vector u(*map);
  u[55] = 1.0;

  // solve
  Epetra_Vector v(*map);
  int ierr = solver->ApplyInverse(u, v);
  CHECK(ierr == 0);

  for (int i = 0; i < 5; i++) CHECK_CLOSE((m->x())[i], v[i], 1e-6);

  delete comm;
};

TEST(SOLVER_FACTORY_WITH_PC) {
  std::cout << "\nChecking solver factory..." << std::endl;

  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(100, 0, *comm));

  Teuchos::ParameterList plist;
  Teuchos::ParameterList& slist = plist.sublist("pcg");
  slist.set<std::string>("iterative method", "pcg");
  slist.set<std::string>("preconditioning method", "identity");

  // create the pcg operator
  Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
  auto solver = AmanziSolvers::createInverse<Matrix,Epetra_Vector,Epetra_Map>(
      "pcg", plist, m);

  // initial guess
  Epetra_Vector u(*map);
  u[55] = 1.0;

  // solve
  Epetra_Vector v(*map);
  int ierr = solver->ApplyInverse(u, v);
  CHECK(ierr == 0);

  for (int i = 0; i < 5; i++) CHECK_CLOSE((m->x())[i], v[i], 1e-6);

  delete comm;
};

TEST(GMRES_WITH_GMRES_PC) {
  std::cout << "\nChecking two-level preconditioner (pcg with gmres pc)..." << std::endl;

  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(100, 0, *comm));

  // create the gmres preconditioner with tolerance 0.01
  Teuchos::Array<std::string> options({ "relative residual" });
  Teuchos::ParameterList plist;
  auto& slist = plist.sublist("gmres");
  slist.set<std::string>("iterative method", "gmres")
       .set<std::string>("preconditioning method", "identity")
       .sublist("gmres parameters").set<double>("error tolerance", 1.0e-3)
                                   .set<Teuchos::Array<std::string>>("convergence criteria", options);

  Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
  auto pc = AmanziSolvers::createInverse<Matrix,Epetra_Vector,Epetra_Map>("gmres", plist, m);

  // create the gmres solver with tolerance 1e-6
  AmanziSolvers::IterativeMethodPCG<Matrix,Amanzi::Matrix<Epetra_Vector,Epetra_Map>,Epetra_Vector,Epetra_Map> solver;
  solver.set_inverse_parameters(plist);
  solver.set_matrices(m, pc);
  solver.set_tolerance(1e-12);
  solver.InitializeInverse();
  solver.ComputeInverse();

  // initial guess
  Epetra_Vector u(*map), v(*map);
  u[55] = 1.0;

  // solve
  int ierr = solver.ApplyInverse(u, v);
  CHECK(ierr == 0);

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
  auto solver = AmanziSolvers::createInverse<Matrix,Epetra_Vector,Epetra_Map>("gmres", plist, m);

  // initial guess
  Epetra_Vector u(*map);
  u[55] = 1.0;

  // solve
  Epetra_Vector v(*map);
  int ierr = solver->ApplyInverse(u, v);
  CHECK(ierr == 0);

  for (int i = 0; i < 5; i++) CHECK_CLOSE((m->x())[i], v[i], 1e-6);

  delete comm;
};

}  // suite SOLVERS




