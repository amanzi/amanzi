/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include "UnitTest++.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "TestReporterStdout.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

#include "errors.hh"

#include "FunctionConstant.hh"
#include "FunctionDistance.hh"
#include "FunctionFactory.hh"
#include "FunctionPolynomial.hh"
#include "FunctionSmoothStep.hh"
#include "FunctionTabular.hh"

using namespace Amanzi;

int
main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize();
  auto status = UnitTest::RunAllTests();
  Kokkos::finalize();
  return status;
}

SUITE(malformed_parameter_list)
{
  TEST(no_function_sublist)
  {
    Teuchos::ParameterList list;
    FunctionFactory fact;
    // Function *f = fact.Create(list);
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
  TEST(unknown_function_sublist)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("fubar");
    FunctionFactory fact;
    // Function *f = fact.Create(list);
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
  TEST(spurious_parameter)
  {
    Teuchos::ParameterList list;
    list.set("fubar", 0);
    Teuchos::ParameterList& sublist = list.sublist("function-constant");
    sublist.set("value", 1.0);
    FunctionFactory fact;
    // Function *f = fact.Create(list);
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
  TEST(extraneous_function_sublist)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist1 = list.sublist("function-constant");
    sublist1.set("value", 1.0);
    // Another sublist.  Note that sublists are processed in alphabetic order.
    Teuchos::ParameterList& sublist2 = list.sublist("zzz-other-function");
    sublist2.set("fubar", 1.0);
    FunctionFactory fact;
    // Function *f = fact.Create(list);
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
}

SUITE(constant_factory)
{
  TEST(create)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-constant");
    sublist.set("value", 2.0);
    FunctionFactory fact;
    Function* f = fact.Create(list);
    Kokkos::View<double*> x("x", 1);
    x(0) = 9.0;
    CHECK_EQUAL((*f)(x), 2.0);
  }
  TEST(missing_parameter)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-constant");
    FunctionFactory fact;
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
}

SUITE(tabular_factory)
{
  TEST(create)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-tabular");
    Teuchos::Array<double> x(2);
    Teuchos::Array<double> y(2);
    x[0] = 0.0;
    x[1] = 1.0;
    y[0] = 2.0;
    y[1] = 3.0;
    sublist.set("x values", x);
    sublist.set("y values", y);
    FunctionFactory fact;
    Function* f = fact.Create(list);
    Kokkos::View<double*> t("t", 1);
    t(0) = 0.5;
    CHECK_EQUAL((*f)(t), 2.5);
  }
  TEST(create_with_row_coordinate)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-tabular");
    Teuchos::Array<double> x(2);
    Teuchos::Array<double> y(2);
    x[0] = 0.0;
    x[1] = 1.0;
    y[0] = 2.0;
    y[1] = 3.0;
    sublist.set("x values", x);
    sublist.set("x coordinate", "t");
    sublist.set("y values", y);
    FunctionFactory fact;
    Function* f = fact.Create(list);
    Kokkos::View<double*> t("t", 1);
    t(0) = 0.5;
    CHECK_EQUAL((*f)(t), 2.5);
  }
  TEST(create_with_form)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-tabular");
    Teuchos::Array<double> x(2);
    Teuchos::Array<double> y(2);
    x[0] = 0.0;
    x[1] = 1.0;
    y[0] = 2.0;
    y[1] = 3.0;
    sublist.set("x values", x);
    sublist.set("y values", y);
    Teuchos::Array<std::string> forms(1, "constant");
    sublist.set("forms", forms);
    FunctionFactory fact;
    Function* f = fact.Create(list);
    Kokkos::View<double*> t("t", 1);
    t(0) = 0.5;
    CHECK_EQUAL(2.0, (*f)(t));
  }
  TEST(missing_parameter)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-tabular");
    Teuchos::Array<double> x(2);
    x[0] = 0.0;
    x[1] = 1.0;
    sublist.set("x values", x);
    FunctionFactory fact;
    // Function *f = fact.Create(list);
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
  TEST(not_enough_points)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-tabular");
    Teuchos::Array<double> x(1);
    Teuchos::Array<double> y(1);
    x[0] = 0.0;
    y[0] = 0.0;
    sublist.set("x values", x);
    sublist.set("y values", y);
    FunctionFactory fact;
    // Function *f = fact.Create(list);
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
  TEST(not_ordered)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-tabular");
    Teuchos::Array<double> x(2);
    Teuchos::Array<double> y(2);
    x[0] = 1.0;
    x[1] = 1.0;
    y[0] = 2.0;
    y[1] = 3.0;
    sublist.set("x values", x);
    sublist.set("y values", y);
    FunctionFactory fact;
    // Function *f = fact.Create(list);
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
  TEST(bad_values)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-tabular");
    Teuchos::Array<double> x(2);
    Teuchos::Array<double> y(3);
    x[0] = 0.0;
    x[1] = 1.0;
    y[0] = 2.0;
    y[1] = 3.0;
    y[2] = 4.0;
    sublist.set("x values", x);
    sublist.set("y values", y);
    FunctionFactory fact;
    // Function *f = fact.Create(list);
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
  TEST(bad_forms)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-tabular");
    Teuchos::Array<double> x(2);
    Teuchos::Array<double> y(2);
    x[0] = 0.0;
    x[1] = 1.0;
    y[0] = 2.0;
    y[1] = 3.0;
    sublist.set("x values", x);
    sublist.set("y values", y);
    Teuchos::Array<std::string> forms(2, "constant");
    sublist.set("forms", forms);
    FunctionFactory fact;
    // Function *f = fact.Create(list);
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
    forms.pop_back();
    forms[0] = "fubar";
    sublist.remove("forms");
    sublist.set("forms", forms);
    // Function *f = fact.Create(list);
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
}

SUITE(bilinear_factory)
{
  TEST(create)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-bilinear");
    sublist.set("file", "test/bilinear.h5");
    sublist.set("row header", "/times");
    sublist.set("row coordinate", "time");
    sublist.set("column header", "/x");
    sublist.set("column coordinate", "x");
    sublist.set("value header", "/values");
    FunctionFactory fact;
    Function* f = fact.Create(list);
    Kokkos::View<double*> t("t", 2);
    t(0) = 2.;
    t(1) = 2.;
    CHECK_EQUAL((*f)(t), 14.);
  }
  TEST(missing_rows)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-bilinear");
    sublist.set("file", "test/bilinear.h5");
    sublist.set("column header", "/x");
    sublist.set("column coordinate", "x");
    sublist.set("value header", "/values");
    FunctionFactory fact;
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
  TEST(missing_columns)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-bilinear");
    sublist.set("file", "test/bilinear.h5");
    sublist.set("row header", "/times");
    sublist.set("row coordinate", "time");
    sublist.set("value header", "/values");
    FunctionFactory fact;
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
  TEST(missing_values)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-bilinear");
    sublist.set("file", "test/bilinear.h5");
    sublist.set("row header", "/times");
    sublist.set("row coordinate", "time");
    sublist.set("column header", "/x");
    sublist.set("column coordinate", "x");
    FunctionFactory fact;
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
  /*
  TEST(not_enough_points)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList &sublist = list.sublist("function-bilinear");
    sublist.set("file", "test/bilinear_missing.h5");  // missing file
    sublist.set("row header", "/times");
    sublist.set("row coordinate", "time");
    sublist.set("column header", "/x");
    sublist.set("column coordinate", "x");
    sublist.set("value header", "/values");
    FunctionFactory fact;
    CHECK_THROW(Function *f = fact.Create(list), Errors::Message);
  }
  TEST(not_ordered)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList &sublist = list.sublist("function-bilinear");
    sublist.set("file", "test/bilinear_unsort.h5");  // missing file
    sublist.set("row header", "/times");
    sublist.set("row coordinate", "time");
    sublist.set("column header", "/x");
    sublist.set("column coordinate", "x");
    sublist.set("value header", "/values");
    FunctionFactory fact;
    CHECK_THROW(Function *f = fact.Create(list), Errors::Message);
  }
  */
}

SUITE(smooth_step_factory)
{
  TEST(create)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-smooth-step");
    sublist.set("x0", 1.0);
    sublist.set("y0", 2.0);
    sublist.set("x1", 3.0);
    sublist.set("y1", 4.0);
    FunctionFactory fact;
    Function* f = fact.Create(list);
    Kokkos::View<double*> x("x", 1);
    x(0) = 2.;
    CHECK_EQUAL((*f)(x), 3.0);
  }
  TEST(missing_parameter)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-smooth-step");
    sublist.set("x0", 1.0);
    sublist.set("y0", 2.0);
    sublist.set("x1", 3.0);
    FunctionFactory fact;
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
  TEST(bad_parameter)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-smooth-step");
    sublist.set("x0", 3.0);
    sublist.set("y0", 2.0);
    sublist.set("x1", 3.0);
    sublist.set("y1", 4.0);
    FunctionFactory fact;
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
}

SUITE(polynomial_factory)
{
  TEST(create)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-polynomial");
    Teuchos::Array<double> c(2);
    Teuchos::Array<int> p(2);
    c[0] = 1.0;
    c[1] = 2.0;
    p[0] = 1;
    p[1] = 0;
    sublist.set("coefficients", c);
    sublist.set("exponents", p);
    FunctionFactory fact;
    Function* f = fact.Create(list);
    Kokkos::View<double*> t("t", 1);
    t(0) = 0.5;
    CHECK_EQUAL((*f)(t), 2.5);
    delete f;
    // Now add the optional reference point argument
    sublist.set("reference point", -1.0);
    f = fact.Create(list);
    t[0] = -0.5;
    CHECK_EQUAL((*f)(t), 2.5);
    delete f;
  }
  TEST(missing_parameter)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-polynomial");
    Teuchos::Array<double> c(2);
    c[0] = 0.0;
    c[1] = 1.0;
    sublist.set("coefficients", c);
    FunctionFactory fact;
    // Function *f = fact.Create(list);
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
  TEST(zero_sized_vectors)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-polynomial");
    Teuchos::Array<double> c(0);
    Teuchos::Array<int> p(0);
    sublist.set("coefficients", c);
    sublist.set("exponents", p);
    FunctionFactory fact;
    // Function *f = fact.Create(list);
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
  TEST(bad_values)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-polynomial");
    Teuchos::Array<double> c(2);
    Teuchos::Array<int> p(3);
    c[0] = 1.0;
    c[1] = 1.0;
    p[0] = 2;
    p[1] = 3;
    p[2] = 4;
    sublist.set("coefficients", c);
    sublist.set("exponents", p);
    FunctionFactory fact;
    // Function *f = fact.Create(list);
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
}

SUITE(linear_factory)
{
  TEST(create)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-linear");
    Teuchos::Array<double> grad(2);
    grad[0] = 1.0;
    grad[1] = 2.0;
    sublist.set("y0", 1.0);
    sublist.set("gradient", grad);
    FunctionFactory fact;
    Function* f = fact.Create(list);
    Kokkos::View<double*> x("x", 2);
    x(0) = 1.;
    x(1) = 2.;
    CHECK_EQUAL(6.0, (*f)(x));
    delete f;
    // Now add the optional x0 parameter
    Teuchos::Array<double> x0(2);
    x0[0] = -1.0;
    x0[1] = 1.0;
    sublist.set("x0", x0);
    f = fact.Create(list);
    CHECK_EQUAL(5.0, (*f)(x));
    delete f;
  }
  TEST(missing_parameter)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-linear");
    sublist.set("y0", 1.0);
    FunctionFactory fact;
    // Function *f = fact.Create(list);
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
    sublist.set("gradient", 1.0);
    // Function *f = fact.Create(list);
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
  TEST(bad_values)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-linear");
    Teuchos::Array<double> grad(2);
    grad[0] = 1.0;
    grad[1] = 2.0;
    sublist.set("y0", 1.0);
    sublist.set("gradient", grad);
    Teuchos::Array<double> x0(3, 1.0);
    sublist.set("x0", x0);
    FunctionFactory fact;
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
}

SUITE(separable_factory)
{
  TEST(create)
  {
    // The parameter list we are creating:
    // <ParameterList name="">
    //   <ParameterList name="function-separable">
    //     <ParameterList name="function1">
    //       <ParameterList name="function-constant">
    //         <Parameter name="value" type="double" value="2.0" />
    //       </ParameterList>
    //     </ParameterList>
    //     <ParameterList name="function2">
    //       <ParameterList name="function-constant">
    //         <Parameter name="value" type="double" value="3.0" />
    //       </ParameterList>
    //     </ParameterList>
    //   </ParameterList>
    // </ParameterList>
    Teuchos::ParameterList list;
    Teuchos::ParameterList& flist = list.sublist("function-separable");
    Teuchos::ParameterList& list1 = flist.sublist("function1");
    Teuchos::ParameterList& flist1 = list1.sublist("function-constant");
    flist1.set("value", 2.0);
    Teuchos::ParameterList& list2 = flist.sublist("function2");
    Teuchos::ParameterList& flist2 = list2.sublist("function-linear");
    Teuchos::Array<double> grad(2);
    grad[0] = 1.0;
    grad[1] = 2.0;
    flist2.set("y0", 0.0);
    flist2.set("gradient", grad);
    FunctionFactory factory;
    Function* f = factory.Create(list);
    Kokkos::View<double*> x("x", 3);
    x(0) = 0.;
    x(1) = 1.;
    x(2) = -1.;
    CHECK_EQUAL(-2.0, (*f)(x));
    delete f;
  }
  TEST(create_nested)
  {
    // The parameter list we are creating:
    // <ParameterList name="">
    //   <ParameterList name="function-separable">
    //     <ParameterList name="function1">
    //       <ParameterList name="function-constant">
    //         <Parameter name="value" type="double" value="2.0" />
    //       </ParameterList>
    //     </ParameterList>
    //     <ParameterList name="function2">
    //       <ParameterList name="function-separable">
    //         <ParameterList name="function1">
    //           <ParameterList name="function-constant">
    //             <Parameter name="value" type="double" value="3.0" />
    //           </ParameterList>
    //         </ParameterList>
    //         <ParameterList name="function2">
    //           <ParameterList name="function-constant">
    //             <Parameter name="value" type="double" value="4.0" />
    //           </ParameterList>
    //         </ParameterList>
    //       </ParameterList>
    //     </ParameterList>
    //   </ParameterList>
    // </ParameterList>
    Teuchos::ParameterList list;
    Teuchos::ParameterList& flist = list.sublist("function-separable");
    Teuchos::ParameterList& listx = flist.sublist("function1");
    Teuchos::ParameterList& flistx = listx.sublist("function-constant");
    flistx.set("value", 2.0);
    Teuchos::ParameterList& listyz = flist.sublist("function2");
    Teuchos::ParameterList& flistyz = listyz.sublist("function-separable");
    Teuchos::ParameterList& listy = flistyz.sublist("function1");
    Teuchos::ParameterList& flisty = listy.sublist("function-constant");
    flisty.set("value", 3.0);
    Teuchos::ParameterList& listz = flistyz.sublist("function2");
    Teuchos::ParameterList& flistz = listz.sublist("function-constant");
    flistz.set("value", 4.0);
    FunctionFactory factory;
    Function* f = factory.Create(list);
    Kokkos::View<double*> x("x", 3);
    x(0) = 0.;
    x(1) = 0.;
    x(2) = 0.;
    CHECK_EQUAL(24.0, (*f)(x));
    delete f;
  }
  TEST(missing_sublist)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& flist = list.sublist("function-separable");
    FunctionFactory factory;
    // Function *f = factory.Create(list);
    CHECK_THROW(Function* f = factory.Create(list), Errors::Message);
    // Now define the function1 sublist.
    Teuchos::ParameterList& list1 = flist.sublist("function1");
    Teuchos::ParameterList& flist1 = list1.sublist("function-constant");
    flist1.set("value", 2.0);
    // Function *f = factory.Create(list);
    CHECK_THROW(Function* f = factory.Create(list), Errors::Message);
  }
}

SUITE(static_head_factory)
{
  TEST(create)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-static-head");
    sublist.set("p0", 1.0);
    sublist.set("density", 4.0);
    sublist.set("gravity", 0.5);
    sublist.set("space dimension", 3);
    sublist.sublist("water table elevation")
      .sublist("function-constant")
      .set("value", 3.0);
    FunctionFactory factory;
    Function* f = factory.Create(list);
    Kokkos::View<double*> x("x", 4);
    x(0) = 0.;
    x(1) = 0.;
    x(2) = 0.;
    x(3) = 1.;
    CHECK_EQUAL(5.0, (*f)(x));
    Kokkos::View<double*> y("y", 4);
    y(0) = 1.;
    y(1) = 1.;
    y(2) = 1.;
    y(3) = 4.;
    CHECK_EQUAL(-1.0, (*f)(y));
  }
  TEST(missing_parameter)
  {
    Teuchos::ParameterList list;
    FunctionFactory factory;
    // Function *f = factory.Create(list);
    CHECK_THROW(Function* f = factory.Create(list), Errors::Message);
    Teuchos::ParameterList& sublist = list.sublist("function-static-head");
    // sublist.set("p0", 1.0);
    sublist.set("density", 4.0).set("gravity", 0.5);
    sublist.sublist("water table elevation")
      .sublist("function-constant")
      .set("value", 3.0);
    // Function *f = factory.Create(list);
    CHECK_THROW(Function* f = factory.Create(list), Errors::Message);
    sublist.set("p0", 1.0);
    sublist.remove("density");
    // Function *f = factory.Create(list);
    CHECK_THROW(Function* f = factory.Create(list), Errors::Message);
    sublist.set("density", 4.0);
    sublist.remove("gravity");
    // Function *f = factory.Create(list);
    CHECK_THROW(Function* f = factory.Create(list), Errors::Message);
    sublist.set("gravity", 0.5);
    sublist.remove("water table elevation");
    // Function *f = factory.Create(list);
    CHECK_THROW(Function* f = factory.Create(list), Errors::Message);
    sublist.sublist("water table elevation")
      .sublist("function-fubar"); // bad function plist
    // Function *f = factory.Create(list);
    CHECK_THROW(Function* f = factory.Create(list), Errors::Message);
  }
}

SUITE(distance_factory)
{
  TEST(create)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-squaredistance");
    Teuchos::Array<double> x0(2), metric(2);
    x0[0] = 1.0;
    x0[1] = 0.0;
    metric[0] = 1.0;
    metric[1] = 1.0;
    sublist.set<Teuchos::Array<double>>("x0", x0).set<Teuchos::Array<double>>(
      "metric", metric);
    FunctionFactory fact;
    Function* f = fact.Create(list);
    Kokkos::View<double*> x("x", 2);
    x(0) = 2.;
    x(1) = 2.;
    CHECK_EQUAL((*f)(x), 5.0);
  }
  TEST(missing_parameter)
  {
    Teuchos::ParameterList list;
    Teuchos::ParameterList& sublist = list.sublist("function-double");
    FunctionFactory fact;
    CHECK_THROW(Function* f = fact.Create(list), Errors::Message);
  }
}
