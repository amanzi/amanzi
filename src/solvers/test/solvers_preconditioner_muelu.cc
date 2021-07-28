/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "MueLu.hpp"
#include "MueLu_TpetraOperator.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"

SUITE(SOLVERS)
{
  using CrsMatrix_type = Tpetra::CrsMatrix<double, int, int>;
  typedef Teuchos::Comm<int> Comm_type;
  using Vector = Tpetra::Vector<double,int,int>;
  using moperator = MueLu::TpetraOperator<double,int,int>; 
  using toperator = Tpetra::Operator<double,int,int>;

  TEST(PRECONDITIONER_MUELU)
  {
    Tpetra::MatrixMarket::Reader<CrsMatrix_type> rd;
    Teuchos::RCP<const Comm_type> comm = Teuchos::DefaultComm<int>::getComm();
    const Teuchos::RCP<CrsMatrix_type> A = rd.readSparseFile("matrix_a.mm",comm); 
    Vector x = Vector(A->getDomainMap());
    Vector y = Vector(A->getRangeMap());
    y.putScalar(0); 
    x.randomize();
    A->apply(x,y);

    Vector x_approx = Vector(x, Teuchos::Copy);
    Vector x_err = Vector(x, Teuchos::Copy);
    Vector x_itr = Vector(x, Teuchos::Copy);
    Vector y_c = Vector(y, Teuchos::Copy);
    Vector y_itr = Vector(y, Teuchos::Copy);
    x_approx.putScalar(0);
    x_itr.putScalar(0);
    std::string xml_plist = "preconditioner_muelu_plist.xml";
    Teuchos::RCP<toperator> op_A = A; 
    Teuchos::RCP<moperator> pc = MueLu::CreateTpetraPreconditioner(op_A,xml_plist); 

    for (int i=0; i!=10; ++i) {
      pc->apply(y_itr, x_itr);
      x_approx.update(1, x_itr, 1);
      x_err.update(1, x, -1, x_approx, 0);
      std::cout<<"x_error/x: "<<x_err.norm2() / x.norm2() <<std::endl;

      A->apply(x_itr, y_c);
      y_itr.update(-1, y_c, 1);
    }
  };

} // suite SOLVERS
