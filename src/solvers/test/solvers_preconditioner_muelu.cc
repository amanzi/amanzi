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

    Vector x_c = Vector(x, Teuchos::Copy);
    x_c.putScalar(0);
    std::string xml_plist = "preconditioner_muelu_plist.xml";
    Teuchos::RCP<toperator> op_A = A; 
    Teuchos::RCP<moperator> pc = MueLu::CreateTpetraPreconditioner(op_A,xml_plist); 
    pc->apply(y, x_c);
    x_c.update(-1,x,1); 
    std::cout<<"x_c/x: "<<x_c.norm2() / x.norm2() <<std::endl;
  };

} // suite SOLVERS
