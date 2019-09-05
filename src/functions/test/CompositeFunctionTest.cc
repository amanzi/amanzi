#include <iostream>

#include "UnitTest++.h"
#include "TestReporterStdout.h"

#include "FunctionConstant.hh"
#include "MultiFunction.hh"
#include "errors.hh"

using namespace Amanzi;

int main (int argc, char *argv[]) {
    Kokkos::initialize(argc,argv); 
    auto status =  UnitTest::RunAllTests ();
    Kokkos::finalize(); 
    return status; 
}

TEST(basic_test) {
  Teuchos::RCP<const Function> f1 = Teuchos::rcp(new FunctionConstant(1.0));
  Teuchos::RCP<const Function> f2 = Teuchos::rcp(new FunctionConstant(2.0));

  std::vector<Teuchos::RCP<const Function> > functions;
  functions.push_back(f1);
  functions.push_back(f2);

  Teuchos::RCP<MultiFunction> fc = Teuchos::rcp(new MultiFunction(functions));
  Kokkos::View<double*> x("x",1); 
  x(0) = 3.0; 
  CHECK_EQUAL((*fc)(x)[0], 1.0);
  CHECK_EQUAL((*fc)(x)[1], 2.0);
}

TEST(factory_test) {
  Teuchos::ParameterList list;

  list.set("number of dofs", 2);

  Teuchos::ParameterList &sublist1 = list.sublist("dof 1 function");
  Teuchos::ParameterList &sublist1A = sublist1.sublist("function-constant");
  sublist1A.set("value", 1.0);

  Teuchos::ParameterList &sublist2 = list.sublist("dof 2 function");
  Teuchos::ParameterList &sublist2A = sublist2.sublist("function-constant");
  sublist2A.set("value", 2.0);

  Teuchos::RCP<MultiFunction> f = Teuchos::rcp(new MultiFunction(list));

  Kokkos::View<double*> x("x",1); 
  x(0) = 3.0; 
  CHECK_EQUAL((*f)(x)[0], 1.0);
  CHECK_EQUAL((*f)(x)[1], 2.0);
}
