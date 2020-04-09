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

#define HAVE_TPETRA_PRECONDITIONERS

#include "UnitTest++.h"

#include "exceptions.hh"
#include "Teuchos_RCP.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "AmanziComm.hh"
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
    using CrsMatrix_type = Tpetra::CrsMatrix<double, int, int>;

   public:
    Matrix() {}
    Matrix(const Map_ptr_type& map) : map_(map){};
    ~Matrix(){};
    Matrix(const Matrix& other) : map_(other.map_) {}

    void Init(const std::string& name)
    {
      Teuchos::ParameterList plist;
      plist.set<std::string>("preconditioner type", name);
      std::string params(name);
      Teuchos::ParameterList& tmp = plist.sublist(params.append(" parameters"));

      // create and assemble the matrix
      A_ = Teuchos::rcp(new CrsMatrix_type(map_, map_, 3));
      for (int i = 0; i < N; i++) {
        int indices[3];
        double values[3] = { double(-i), double(2 * i + 1), double(-i - 1) };
        for (int k = 0; k < 3; k++) indices[k] = i + k - 1;
        A_->insertLocalValues(i, 3, values, indices);
      }
      A_->fillComplete(map_, map_);

      // create and assemble the preconditioner
      if (name == "diagonal") {
        preconditioner_ =
          Teuchos::rcp(new PreconditionerDiagonal<Matrix_type, Vector_type>());
      } else if (name == "identity") {
        preconditioner_ =
          Teuchos::rcp(new PreconditionerIdentity<Matrix_type, Vector_type>());
      } else if (name == "ml") {
        PreconditionerFactory<Matrix_type, Vector_type> factory;
        tmp.set<int>("coarse: max size", 5);
        tmp.set<int>("cycle applications", 2);
        tmp.set<int>("ML output", 0);
        preconditioner_ = factory.Create(plist);
      } else {
        tmp.set<int>("max coarse size", 5);
        tmp.set<int>("cycle applications", 1);
        tmp.set<int>("verbosity", 0);
        PreconditionerFactory<Matrix_type, Vector_type> factory;
        preconditioner_ = factory.Create(plist);
      }
      preconditioner_->Update(A_);
    };

    int apply(const Vector_type& v, Vector_type& mv) const
    {
      //    std::cout << "In Matrix::Apply: v = " << v.norm2();
      A_->apply(v, mv);
      //    std::cout << " mv = " << mv.norm2() << std::endl;
      return 0;
    }
    int applyInverse(const Vector_type& v, Vector_type& hv) const
    {
      return preconditioner_->applyInverse(v, hv);
    }

    const Map_ptr_type& getDomainMap() const { return map_; }
    const Map_ptr_type& getRangeMap() const { return map_; }

   private:
    Map_ptr_type map_;
    Teuchos::RCP<CrsMatrix_type> A_;
    Teuchos::RCP<Preconditioner<Matrix_type, Vector_type>> preconditioner_;
  };


  TEST(DIAGONAL_PRECONDITIONER)
  {
    std::cout << "\nComparison of preconditioners for N=125" << std::endl;

    auto comm = getDefaultComm();
    auto map = Teuchos::rcp(new Map_type(N, 0, comm));

    // create the pcg operator
    Teuchos::RCP<Matrix> m = Teuchos::rcp(new Matrix(map));
    AmanziSolvers::LinearOperatorPCG<Matrix, Vector_type, Map_type> pcg(m, m);
    pcg.Init();
    pcg.set_tolerance(1e-12);
    pcg.set_max_itrs(200);

    Vector_type u(map), v(map);
    {
      auto uv = u.getLocalViewHost();
      for (int i = 0; i < N; i++) uv(i, 0) = 1.0 / (i + 2.0);
    }
    u.sync_device();

    // solving with preconditioner
    std::vector<std::string> prec_names{ "identity", "diagonal" };
    for (auto name : prec_names) {
      m->Init(name);
      v.putScalar(0.0);
      std::cout << "Preconditioner: " << name << std::endl;
      pcg.applyInverse(u, v);

      v.sync_host();
      {
        auto vv = v.getLocalViewHost();
        CHECK_CLOSE(11.03249773994628, vv(0, 0), 1e-6);
        CHECK_CLOSE(10.53249773994628, vv(1, 0), 1e-6);
      }
    }
  };
}
