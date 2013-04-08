#include "UnitTest++.h"
#include "TestReporterStdout.h"

#include "constant-function.hh"
#include "vector_function_factory.hh"
#include "vector_function.hh"
#include "composite_function.hh"
#include "errors.hh"

#include <iostream>

using namespace Amanzi;

int main (int argc, char *argv[]) {
    return UnitTest::RunAllTests ();
}

TEST(basic_test) {
  Teuchos::RCP<const Function> f1 = Teuchos::rcp(new ConstantFunction(1.0));
  Teuchos::RCP<const Function> f2 = Teuchos::rcp(new ConstantFunction(2.0));

  std::vector<Teuchos::RCP<const Function> > functions;
  functions.push_back(f1);
  functions.push_back(f2);

  Teuchos::RCP<VectorFunction> fc = Teuchos::rcp(new CompositeFunction(functions));
  double x = 3.0;
  CHECK_EQUAL((*fc)(&x)[0], 1.0);
  CHECK_EQUAL((*fc)(&x)[1], 2.0);
}

TEST(factory_test) {
  Teuchos::ParameterList list;

  list.set("Number of DoFs", 2);
  list.set("Function type", "composite function");

  Teuchos::ParameterList &sublist1 = list.sublist("DoF 1 Function");
  Teuchos::ParameterList &sublist1A = sublist1.sublist("function-constant");
  sublist1A.set("value", 1.0);

  Teuchos::ParameterList &sublist2 = list.sublist("DoF 2 Function");
  Teuchos::ParameterList &sublist2A = sublist2.sublist("function-constant");
  sublist2A.set("value", 2.0);

  VectorFunctionFactory factory;
  Teuchos::RCP<VectorFunction> f = factory.Create(list);

  double x = 3.0;
  CHECK_EQUAL((*f)(&x)[0], 1.0);
  CHECK_EQUAL((*f)(&x)[1], 2.0);
}
