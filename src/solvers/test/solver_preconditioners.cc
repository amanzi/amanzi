#include <iostream>
#include <string>

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "exceptions.hh"
#include "LinearOperatorPCG.hh"
#include "PreconditionerDiagonal.hh"
#include "PreconditionerIdentity.hh"

SUITE(SOLVERS) {
const int N = 25;
using namespace Amanzi;
using namespace Amanzi::AmanziPreconditioners;

class Matrix {
 public:
  Matrix() {}
  Matrix(const Teuchos::RCP<Epetra_Map>& map) : map_(map) {};
  ~Matrix() {};
  Matrix(const Matrix& other) : map_(other.map_) {}

  void Init(std::string& name, const Epetra_Map& map) {
    if (name == "diagonal") {
      preconditioner_ = Teuchos::rcp(new PreconditionerDiagonal());
    }
    A_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map, map, 3));
    for (int i = 0; i < N; i++) {
      int indices[3];
      double values[3] = {-i, 2 * i + 1, -i - 1};
      for (int k = 0; k < 3; k++) indices[k] = i + k - 1; 
      A_->InsertMyValues(i, 3, values, indices);
    }
    A_->FillComplete(map, map);
    preconditioner_->Update(A_);
  };    

  virtual int Apply(const Epetra_Vector& v, Epetra_Vector& mv) const { 
    A_->Apply(v, mv);
  }
  virtual int ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const {
    preconditioner_->ApplyInverse(v, hv);
  }

  virtual const Epetra_Map& DomainMap() const { return *map_; }
  virtual const Epetra_Map& RangeMap() const { return *map_; }

 private:
  Teuchos::RCP<Epetra_Map> map_;
  Teuchos::RCP<Epetra_CrsMatrix> A_;
  Teuchos::RCP<Preconditioner> preconditioner_;
};

TEST(DIAGONAL_PRECONDITONER) {
  std::cout << "Checking PCG solver..." << std::endl;

  Epetra_MpiComm* comm = new Epetra_MpiComm(MPI_COMM_SELF);
  Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(N, 0, *comm));

  // create the pcg operator
  Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
  std::string prec_name("diagonal");
  m->Init(prec_name, *map);

  AmanziSolvers::LinearOperatorPCG<Matrix, Epetra_Vector, Epetra_Map> pcg(m, m);
  pcg.Init();

  // initial guess
  Epetra_Vector u(*map);
  u[0] = 1.0;

  // solve
  Epetra_Vector v(*map);
  pcg.ApplyInverse(u, v);

  CHECK_CLOSE(3.81595818e+0, v[0], 1e-6);
  CHECK_CLOSE(2.81595818e+0, v[1], 1e-6);

  delete comm;
};

}




