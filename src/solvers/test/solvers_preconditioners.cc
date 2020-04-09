/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <iostream>
#include <string>

#ifdef _OPENMP
#  include "omp.h"
#endif

#define HAVE_EPETRA_PRECONDITIONERS
#define HAVE_TRILINOS_PRECONDITIONERS
#define HAVE_HYPRE_PRECONDITIONERS

#include "Teuchos_RCP.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "UnitTest++.h"

#include "exceptions.hh"
#include "LinearOperatorPCG.hh"
#include "PreconditionerFactory.hh"
#include "PreconditionerDiagonal.hh"
#include "PreconditionerIdentity.hh"

SUITE(SOLVERS)
{
  const int N = 125;
  using namespace Amanzi;
  using namespace Amanzi::AmanziPreconditioners;

  class Matrix {
   public:
    Matrix() {}
    Matrix(const Teuchos::RCP<Epetra_Map>& map) : map_(map){};
    ~Matrix(){};
    Matrix(const Matrix& other) : map_(other.map_) {}

    void Init(const std::string& name, const Epetra_Map& map)
    {
      Teuchos::ParameterList plist;
      plist.set<std::string>("preconditioner type", name);
      std::string params(name);
      Teuchos::ParameterList& tmp = plist.sublist(params.append(" parameters"));

      if (name == "diagonal") {
        preconditioner_ = Teuchos::rcp(
          new PreconditionerDiagonal<Epetra_RowMatrix, Epetra_MultiVector>());
      } else if (name == "identity") {
        preconditioner_ = Teuchos::rcp(
          new PreconditionerIdentity<Epetra_RowMatrix, Epetra_MultiVector>());
      } else if (name == "ml") {
        PreconditionerFactory<Epetra_RowMatrix, Epetra_MultiVector> factory;
        tmp.set<int>("coarse: max size", 5);
        tmp.set<int>("cycle applications", 2);
        tmp.set<int>("ML output", 0);
        preconditioner_ = factory.Create(plist);
      } else {
        tmp.set<int>("max coarse size", 5);
        tmp.set<int>("cycle applications", 1);
        tmp.set<int>("verbosity", 0);
        PreconditionerFactory<Epetra_RowMatrix, Epetra_MultiVector> factory;
        preconditioner_ = factory.Create(plist);
      }
      A_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map, map, 3));
      for (int i = 0; i < N; i++) {
        int indices[3];
        double values[3] = { double(-i), double(2 * i + 1), double(-i - 1) };
        for (int k = 0; k < 3; k++) indices[k] = i + k - 1;
        A_->InsertMyValues(i, 3, values, indices);
      }
      A_->FillComplete(map, map);
      preconditioner_->update(A_);
    };

    virtual int Apply(const Epetra_MultiVector& v, Epetra_MultiVector& mv) const
    {
      return A_->apply(v, mv);
    }
    virtual int
    ApplyInverse(const Epetra_MultiVector& v, Epetra_MultiVector& hv) const
    {
      return preconditioner_->applyInverse(v, hv);
    }

    virtual const Epetra_Map& DomainMap() const { return *map_; }
    virtual const Epetra_Map& RangeMap() const { return *map_; }

   private:
    Teuchos::RCP<Epetra_Map> map_;
    Teuchos::RCP<Epetra_CrsMatrix> A_;
    Teuchos::RCP<Preconditioner<Epetra_RowMatrix, Epetra_MultiVector>>
      preconditioner_;
  };


  TEST(DIAGONAL_PRECONDITIONER)
  {
    std::cout << "\nComparison of preconditioners for N=125" << std::endl;

    Epetra_MpiComm comm(MPI_COMM_SELF);
    Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(N, 0, comm));

    // create the pcg operator
    Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
    AmanziSolvers::LinearOperatorPCG<Matrix, Epetra_Vector, Epetra_Map> pcg(m,
                                                                            m);
    pcg.Init();
    pcg.set_tolerance(1e-12);
    pcg.set_max_itrs(200);

    Epetra_Vector u(*map), v(*map);
    for (int i = 0; i < N; i++) u[i] = 1.0 / (i + 2.0);

    // solving with identity preconditioner
    std::string prec_names[4];
    prec_names[0] = "identity";
    prec_names[1] = "diagonal";
    prec_names[2] = "boomer amg";
    prec_names[3] = "ml";
    for (int n = 0; n < 4; n++) {
      m->Init(prec_names[n], *map);

      v.putScalar(0.0);
      printf("Preconditioner: %s\n", prec_names[n].c_str());
      pcg.applyInverse(u, v);

      CHECK_CLOSE(11.03249773994628, v[0], 1e-6);
      CHECK_CLOSE(10.53249773994628, v[1], 1e-6);
    }
  };


#ifdef _OPENMP
  TEST(DIAGONAL_PRECONDITIONER_OPENMP)
  {
    std::cout
      << "\nComparison of preconditioners for N=125 using OpenMP directives"
      << std::endl;

    Epetra_MpiComm comm(MPI_COMM_SELF);
    Teuchos::RCP<Epetra_Map> map = Teuchos::rcp(new Epetra_Map(N, 0, comm));

    Teuchos::ParameterList plist;
    plist.sublist("verbose object").set<std::string>("verbosity level", "low");

    // solving with identity preconditioner
    std::string prec_names[4];
    prec_names[0] = "identity";
    prec_names[1] = "diagonal";

    // parallel run
    double cpu0 = omp_get_wtime();

#  pragma omp parallel for shared(map) num_threads(2)
    for (int n = 0; n < 2; n++) {
      Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
      m->Init(prec_names[n], *map);

      AmanziSolvers::LinearOperatorPCG<Matrix, Epetra_Vector, Epetra_Map> pcg(
        m, m);
      pcg.Init(plist);
      pcg.set_tolerance(1e-12);
      pcg.set_max_itrs(200);

      Epetra_Vector u(*map), v(*map);
      for (int i = 0; i < N; i++) u[i] = 1.0 / (i + 2.0);

      pcg.applyInverse(u, v);
    }

    int nthreads = omp_get_max_threads();
    double cpu1 = omp_get_wtime();
    std::cout << "CPU (parallel): " << cpu1 - cpu0
              << " [sec]  threads=" << nthreads << std::endl;

    // serial run
    cpu0 = omp_get_wtime();

    for (int n = 0; n < 2; n++) {
      Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
      m->Init(prec_names[n], *map);

      AmanziSolvers::LinearOperatorPCG<Matrix, Epetra_Vector, Epetra_Map> pcg(
        m, m);
      pcg.Init(plist);
      pcg.set_tolerance(1e-12);
      pcg.set_max_itrs(200);

      Epetra_Vector u(*map), v(*map);
      for (int i = 0; i < N; i++) u[i] = 1.0 / (i + 2.0);

      pcg.applyInverse(u, v);
    }

    cpu1 = omp_get_wtime();
    std::cout << "CPU (serial):   " << cpu1 - cpu0 << " [sec] " << std::endl;
  };
#endif
}
