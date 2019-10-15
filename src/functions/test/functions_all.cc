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
#include "TestReporterStdout.h"

#include "FunctionConstant.hh"
#include "FunctionSmoothStep.hh"
#include "FunctionTabular.hh"
#include "FunctionPolynomial.hh"
#include "FunctionLinear.hh"
#include "FunctionSeparable.hh"
#include "FunctionStaticHead.hh"
#include "FunctionBilinear.hh"
#include "errors.hh"
#include "HDF5Reader.hh"

using namespace Amanzi;

int
main(int argc, char* argv[])
{
  Kokkos::initialize();
  int status = UnitTest::RunAllTests();
  Kokkos::finalize();
  return status;
}

TEST(constant_test)
{
  Function* f = new FunctionConstant(1.0);
  Kokkos::View<double*> x("x", 2);
  x(0) = 1;
  x(1) = 3.0;
  CHECK_EQUAL((*f)(x), 1.0);
  Function* g = f->Clone();
  delete f;
  CHECK_EQUAL((*g)(x), 1.0);
  delete g;
}

TEST(smooth_step_test)
{
  Function* f = new FunctionSmoothStep(1.0, 1.0, 3.0, 0.0);
  Kokkos::View<double*> x("x", 2);
  x(0) = 1;
  x(1) = 0.0;
  x[0] = 0.0;
  CHECK_EQUAL((*f)(x), 1.0);
  x[0] = 1.0;
  CHECK_EQUAL((*f)(x), 1.0);
  x[0] = 1.5;
  CHECK_EQUAL((*f)(x), 0.84375);
  x[0] = 2.0;
  CHECK_EQUAL((*f)(x), 0.5);
  x[0] = 2.5;
  CHECK_EQUAL((*f)(x), 0.15625);
  x[0] = 3.0;
  CHECK_EQUAL((*f)(x), 0.0);
  x[0] = 4.0;
  CHECK_EQUAL((*f)(x), 0.0);
  Function* g = f->Clone();
  delete f;
  x[0] = 2.0;
  CHECK_EQUAL((*g)(x), 0.5);
  delete g;
  CHECK_THROW(FunctionSmoothStep f(3.0, 1.0, 3.0, 0.0), Errors::Message);
}


TEST(tabular_test)
{
  Kokkos::View<double*> x("x", 5);
  x(0) = 0.0;
  x(1) = 1.0;
  x(2) = 3.0;
  x(3) = 3.5;
  x(4) = 3.5;
  Kokkos::View<double*> y("y", 5);
  y(0) = 1.0;
  y(1) = 3.0;
  y(2) = 2.0;
  y(3) = 0.0;
  y(4) = 0.0;
  Kokkos::View<double*> xvec = Kokkos::subview(x, Kokkos::make_pair(0, 4));
  Kokkos::View<double*> yvec = Kokkos::subview(y, Kokkos::make_pair(0, 4));
  int xi = 0;
  Function* f = new FunctionTabular(xvec, yvec, xi);
  Kokkos::View<double*> z("z", 1);
  z(0) = -1;
  CHECK_EQUAL((*f)(z), 1.0);
  z[0] = 0.0;
  CHECK_EQUAL((*f)(z), 1.0);
  z[0] = 0.5;
  CHECK_EQUAL((*f)(z), 2.0);
  z[0] = 1.0;
  CHECK_EQUAL((*f)(z), 3.0);
  z[0] = 2.0;
  CHECK_EQUAL((*f)(z), 2.5);
  z[0] = 3.0;
  CHECK_EQUAL((*f)(z), 2.0);
  z[0] = 3.25;
  CHECK_EQUAL((*f)(z), 1.0);
  z[0] = 3.5;
  CHECK_EQUAL((*f)(z), 0.0);
  z[0] = 4.0;
  CHECK_EQUAL((*f)(z), 0.0);
  Function* g = f->Clone();
  delete f;
  z[0] = 3.0;
  CHECK_EQUAL((*g)(z), 2.0);
  delete g;

  // Now try with optional form argument
  Kokkos::View<FunctionTabular::Form*> form("form", 3);
  form(0) = FunctionTabular::CONSTANT;
  form(1) = FunctionTabular::LINEAR;
  form(2) = FunctionTabular::CONSTANT;
  Kokkos::View<FunctionTabular::Form*> formvec =
    Kokkos::subview(form, Kokkos::make_pair(0, 3));

  f = new FunctionTabular(xvec, yvec, xi, formvec);
  z[0] = -1.0;
  CHECK_EQUAL((*f)(z), 1.0);
  z[0] = 0.0;
  CHECK_EQUAL((*f)(z), 1.0);
  z[0] = 0.5;
  CHECK_EQUAL((*f)(z), 1.0);
  z[0] = 1.0;
  CHECK_EQUAL((*f)(z), 1.0);
  z[0] = 2.0;
  CHECK_EQUAL((*f)(z), 2.5);
  z[0] = 3.0;
  CHECK_EQUAL((*f)(z), 2.0);
  z[0] = 3.25;
  CHECK_EQUAL((*f)(z), 2.0);
  z[0] = 3.5;
  CHECK_EQUAL((*f)(z), 2.0);
  z[0] = 4.0;
  CHECK_EQUAL((*f)(z), 0.0);
  g = f->Clone();
  delete f;
  z[0] = 3.0;
  CHECK_EQUAL((*g)(z), 2.0);
  delete g;
  // Verify the constructor fails with only a single table entry.
  xvec = Kokkos::subview(x, Kokkos::make_pair(0, 1));
  yvec = Kokkos::subview(y, Kokkos::make_pair(0, 1));
  // f = new FunctionTabular(xvec, yvec);
  CHECK_THROW(f = new FunctionTabular(xvec, yvec, xi), Errors::Message);
  // Verify the constructor fails when the x values are not strictly increasing.
  xvec = Kokkos::subview(x, Kokkos::make_pair(0, 5));
  yvec = Kokkos::subview(y, Kokkos::make_pair(0, 5));
  CHECK_THROW(f = new FunctionTabular(xvec, yvec, xi), Errors::Message);
  // Verify the constructor fails with different number of x and y values.
  xvec = Kokkos::subview(x, Kokkos::make_pair(0, 4));
  CHECK_THROW(f = new FunctionTabular(xvec, yvec, xi), Errors::Message);
  // Verify the constructor fails with the wrong number of form values.
  yvec = Kokkos::subview(
    y, Kokkos::make_pair(0, 4)); // same number of x and y values now
  formvec = Kokkos::subview(formvec, Kokkos::make_pair(0, 2));
  CHECK_THROW(f = new FunctionTabular(xvec, yvec, xi, formvec),
              Errors::Message);
}


SUITE(polynomial_test)
{
  TEST(poly1)
  {
    // polynomial with all positive powers
    Kokkos::View<double*> c("c", 2);
    c(0) = -1.0;
    c(1) = 2.0;
    Kokkos::View<int*> p("p", 2);
    p(0) = 3;
    p(1) = 1;
    Kokkos::View<double*> cvec = Kokkos::subview(c, Kokkos::make_pair(0, 2));
    Kokkos::View<int*> pvec = Kokkos::subview(p, Kokkos::make_pair(0, 2));
    Function* f = new FunctionPolynomial(cvec, pvec);
    Kokkos::View<double*> x("x", 1);
    x(0) = 2.0;
    CHECK_EQUAL((*f)(x), -4.0);
    x[0] = -2.0;
    CHECK_EQUAL((*f)(x), 4.0);
    delete f;
    // same polynomial shifted 2 to the right
    f = new FunctionPolynomial(cvec, pvec, 2.0);
    x[0] = 4.0;
    CHECK_EQUAL((*f)(x), -4.0);
    x[0] = 0.0;
    CHECK_EQUAL((*f)(x), 4.0);
    delete f;
  }
  TEST(poly2)
  {
    // polynomial with positive powers, constant term, and repeated power
    Kokkos::View<double*> c("c", 4);
    c(0) = -1.0;
    c(1) = 4.0;
    c(2) = 3.0;
    c(3) = -1.0;
    Kokkos::View<int*> p("p", 4);
    p(0) = 3;
    p(1) = 0;
    p(2) = 1;
    p(3) = 1;
    Kokkos::View<double*> cvec = Kokkos::subview(c, Kokkos::make_pair(0, 4));
    Kokkos::View<int*> pvec = Kokkos::subview(p, Kokkos::make_pair(0, 4));
    Function* f = new FunctionPolynomial(cvec, pvec);
    Kokkos::View<double*> x("x", 1);
    x(0) = 2.0;
    CHECK_EQUAL((*f)(x), 0.0);
    x[0] = 0.0;
    CHECK_EQUAL((*f)(x), 4.0);
    x[0] = -2.0;
    CHECK_EQUAL((*f)(x), 8.0);
    delete f;
  }
  TEST(poly3)
  {
    // polynomial with negative powers only
    // double c[2] = {2.0, -1.0};
    Kokkos::View<double*> c("c", 2);
    c(0) = 2.0;
    c(1) = -1.0;
    // int p[2] = {-3, -2};
    Kokkos::View<int*> p("p", 2);
    p(0) = -3;
    p(1) = -2;
    // std::vector<double> cvec(c, c+2);
    Kokkos::View<double*> cvec = Kokkos::subview(c, Kokkos::make_pair(0, 2));
    // std::vector<int> pvec(p, p+2);
    Kokkos::View<int*> pvec = Kokkos::subview(p, Kokkos::make_pair(0, 2));
    Function* f = new FunctionPolynomial(cvec, pvec, 1.0);
    // std::vector<double> x(1,3.0);
    Kokkos::View<double*> x("x", 1);
    x(0) = 3.0;
    CHECK_EQUAL((*f)(x), 0.0);
    x[0] = 2.0;
    CHECK_EQUAL((*f)(x), 1.0);
    x[0] = -1.0;
    CHECK_EQUAL((*f)(x), -0.5);
    delete f;
  }
  TEST(poly4)
  {
    // polynomial with positive and negative powers
    Kokkos::View<double*> c("c", 4);
    c(0) = 2.0;
    c(1) = -1.0;
    c(2) = -2.0;
    c(3) = 4.0;
    Kokkos::View<int*> p("p", 4);
    p(0) = 1;
    p(1) = 3;
    p(2) = 0;
    p(3) = -2;
    Kokkos::View<double*> cvec = Kokkos::subview(c, Kokkos::make_pair(0, 4));
    Kokkos::View<int*> pvec = Kokkos::subview(p, Kokkos::make_pair(0, 4));
    Function* f = new FunctionPolynomial(cvec, pvec);
    Kokkos::View<double*> x("x", 1);
    x(0) = 2.0;
    CHECK_EQUAL((*f)(x), -5.0);
    x[0] = 1.0;
    CHECK_EQUAL((*f)(x), 3.0);
    x[0] = -1.0;
    CHECK_EQUAL((*f)(x), 1.0);
    x[0] = -2.0;
    CHECK_EQUAL((*f)(x), 3.0);
    delete f;
  }
  TEST(poly_clone)
  {
    // polynomial with positive and negative powers
    Kokkos::View<double*> c("c", 4);
    c(0) = 2.0;
    c(1) = -1.0;
    c(2) = -2.0;
    c(3) = 4.0;
    Kokkos::View<int*> p("p", 4);
    p(0) = 1;
    p(1) = 3;
    p(2) = 0;
    p(3) = -2;
    Kokkos::View<double*> cvec = Kokkos::subview(c, Kokkos::make_pair(0, 4));
    Kokkos::View<int*> pvec = Kokkos::subview(p, Kokkos::make_pair(0, 4));
    FunctionPolynomial* f = new FunctionPolynomial(cvec, pvec);
    Kokkos::View<double*> x("x", 1);
    x(0) = 2.0;
    CHECK_EQUAL((*f)(x), -5.0);
    FunctionPolynomial* g = f->Clone();
    delete f;
    CHECK_EQUAL((*g)(x), -5.0);
    delete g;
  }
}


TEST(linear_test)
{
  double y0 = 1.0;
  Kokkos::View<double*> grad("grad", 3);
  grad(0) = 1.0;
  grad(1) = 2.0;
  grad(2) = -3.0;
  Kokkos::View<double*> gradvec =
    Kokkos::subview(grad, Kokkos::make_pair(0, 3));
  Function* f;
  // Three variable linear function.
  f = new FunctionLinear(y0, gradvec);
  Kokkos::View<double*> a("a", 3);
  a(0) = 3.;
  a(1) = 2.;
  a(2) = 1.;
  Kokkos::View<double*> b("b", 3);
  b(0) = 1.;
  b(1) = -1.;
  b(2) = 1.;
  CHECK_EQUAL(5.0, (*f)(a));
  CHECK_EQUAL(-3.0, (*f)(b));
  Function* g = f->Clone();
  delete f;
  CHECK_EQUAL(-3.0, (*g)(b));
  delete g;
  // Two variable linear function with x0.
  Kokkos::View<double*> x0("x0", 3);
  x0(0) = 1.0;
  x0(1) = 2.0;
  x0(2) = 3.0;
  Kokkos::View<double*> x0vec = Kokkos::subview(x0, Kokkos::make_pair(0, 2));
  gradvec = Kokkos::subview(grad, Kokkos::make_pair(0, 2));
  f = new FunctionLinear(y0, gradvec, x0vec);
  CHECK_EQUAL(3.0, (*f)(a));
  CHECK_EQUAL(-5.0, (*f)(b));
  // -- check error in case of point smaller than gradient
  Kokkos::View<double*> c("c", 1);
  c(0) = 1.;
  CHECK_THROW((*f)(c), Errors::Message);
  delete f;

  // Single variable linear function.
  gradvec = Kokkos::subview(grad, Kokkos::make_pair(0, 1));
  f = new FunctionLinear(y0, gradvec);
  CHECK_EQUAL(4.0, (*f)(a));
  CHECK_EQUAL(2.0, (*f)(b));
  delete f;
  // Check exception handling.
  gradvec = Kokkos::subview(grad, Kokkos::make_pair(0, 0));
  CHECK_THROW(f = new FunctionLinear(y0, gradvec), Errors::Message);
  CHECK_THROW(f = new FunctionLinear(y0, gradvec, x0vec), Errors::Message);
  gradvec = Kokkos::subview(grad, Kokkos::make_pair(0, 1));
  CHECK_THROW(f = new FunctionLinear(y0, gradvec, x0vec), Errors::Message);
}


SUITE(separable_test)
{
  TEST(separable1)
  {
    // First function a smooth step function
    FunctionSmoothStep f1(1.0, 1.0, 3.0, 0.0);
    // Second function a linear function
    Kokkos::View<double*> grad("grad", 3);
    grad(0) = 1.0;
    grad(1) = 2.0;
    grad(2) = 3.0;
    Kokkos::View<double*> gradvec =
      Kokkos::subview(grad, Kokkos::make_pair(0, 3));
    Function* f2 = new FunctionLinear(1.0, gradvec);
    Function* f = new FunctionSeparable(f1, *f2);
    delete f2;
    Kokkos::View<double*> x("x", 4);
    x(0) = 0.;
    x(1) = 1.;
    x(2) = 2.;
    x(3) = -1.;
    CHECK_EQUAL(3.0, (*f)(x));
    x[0] = 5.0;
    CHECK_EQUAL(0.0, (*f)(x));
    delete f;
  }
  TEST(separable2)
  {
    // separable built from separable
    Function* fx = new FunctionSmoothStep(0.0, 1.0, 1.0, 2.0);
    Function* fy = new FunctionSmoothStep(0.0, 1.0, 1.0, 3.0);
    Function* fz = new FunctionSmoothStep(0.0, 1.0, 1.0, 4.0);
    Function* fyz = new FunctionSeparable(*fy, *fz);
    Function* fxyz = new FunctionSeparable(*fx, *fyz);
    delete fx;
    delete fy;
    delete fz;
    delete fyz;
    Kokkos::View<double*> x0("x0", 3);
    x0(0) = 0.;
    x0(1) = 0.;
    x0(2) = 0.;
    CHECK_EQUAL(1.0, (*fxyz)(x0));
    x0[0] = 1.;
    CHECK_EQUAL(2.0, (*fxyz)(x0));
    x0[0] = 0.;
    x0[1] = 1.;
    CHECK_EQUAL(3.0, (*fxyz)(x0));
    x0[1] = 0.;
    x0[2] = 1.;
    CHECK_EQUAL(4.0, (*fxyz)(x0));
    Kokkos::View<double*> x4("x4", 3);
    x4(0) = 0.5;
    x4(1) = 0.5;
    x4(2) = 0.5;
    CHECK_EQUAL(7.5, (*fxyz)(x4));
    delete fxyz;
  }
}

TEST(static_head_test)
{
  // txy-linear height function
  double y0 = 1.0;
  Kokkos::View<double*> grad("grad", 3);
  grad(0) = 0.0;
  grad(1) = 2.0;
  grad(2) = 3.0;
  Kokkos::View<double*> gradvec =
    Kokkos::subview(grad, Kokkos::make_pair(0, 3));
  Function* h = new FunctionLinear(y0, gradvec);
  double patm = -1.0, rho = 0.5, g = 4.0;
  Function* f = new FunctionStaticHead(patm, rho, g, *h, 3);
  Kokkos::View<double*> a("a", 4);
  a(0) = 0.;
  a(1) = 1.;
  a(2) = 0.5;
  a(3) = 2.;
  CHECK_EQUAL(patm + rho * g * (y0 + grad[1] * a[1] + grad[2] * a[2] - a[3]),
              (*f)(a));
  Kokkos::View<double*> b("b", 4);
  b(0) = 1.;
  b(1) = 2.;
  b(2) = 1.;
  b(3) = 1.;
  CHECK_EQUAL(patm + rho * g * (y0 + grad[1] * b[1] + grad[2] * b[2] - b[3]),
              (*f)(b));
  Function* ff = f->Clone();
  double val = (*f)(b);
  delete f;
  CHECK_EQUAL(val, (*ff)(b));
}

// SUITE(separable_test) {
//  TEST(separable1)
//  {
//    // First function a smooth step function
//    std::auto_ptr<Function> f1(new FunctionSmoothStep(1.0, 1.0, 3.0, 0.0));
//    // Second function a linear function
//    double grad[3] = {1.0, 2.0, 3.0};
//    std::vector<double> gradvec(grad, grad+3);
//    std::auto_ptr<Function> f2(new FunctionLinear(1.0, gradvec));
//    Function *f = new FunctionSeparable(f1, f2);
//    double x[4] = {0.0, 1.0, 2.0, -1.0}; // t, x, y, z
//    CHECK_EQUAL(3.0, (*f)(x));
//    x[0] = 5.0;
//    CHECK_EQUAL(0.0, (*f)(x));
//    delete f;
//  }
//
//  TEST(separable2)
//  {
//    // separable built from separable
//    std::auto_ptr<Function> fx(new FunctionSmoothStep(0.0, 1.0, 1.0, 2.0));
//    std::auto_ptr<Function> fy(new FunctionSmoothStep(0.0, 1.0, 1.0, 3.0));
//    std::auto_ptr<Function> fz(new FunctionSmoothStep(0.0, 1.0, 1.0, 4.0));
//    std::auto_ptr<Function> fyz(new FunctionSeparable(fy, fz));
//    Function *fxyz = new FunctionSeparable(fx, fyz);
//    double x0[3] = {0.0, 0.0, 0.0};
//    double x1[3] = {1.0, 0.0, 0.0};
//    double x2[3] = {0.0, 1.0, 0.0};
//    double x3[3] = {0.0, 0.0, 1.0};
//    double x4[3] = {0.5, 0.5, 0.5};
//    CHECK_EQUAL(1.0, (*fxyz)(x0));
//    CHECK_EQUAL(2.0, (*fxyz)(x1));
//    CHECK_EQUAL(3.0, (*fxyz)(x2));
//    CHECK_EQUAL(4.0, (*fxyz)(x3));
//    CHECK_EQUAL(7.5, (*fxyz)(x4));
//    delete fxyz;
//  }
//}

TEST(bilinear_test)
{
  std::string filename = "test/bilinear.h5";
  HDF5Reader reader(filename);
  std::string row_name = "times";
  // std::vector<double> vec_x;
  Kokkos::View<double*> vec_x;
  int xi = 0;
  reader.ReadData(row_name, vec_x);
  std::string col_name = "x";
  // std::vector<double> vec_y;
  Kokkos::View<double*> vec_y;
  int yi = 1;
  reader.ReadData(col_name, vec_y);
  std::string v_name = "values";

  Kokkos::View<double**> mat_v;
  reader.ReadMatData(v_name, mat_v);
  Function* f = new FunctionBilinear(vec_x, vec_y, mat_v, xi, yi);
  // Corners
  // std::vector<double> z(2,0.);
  Kokkos::View<double*> z("z", 2);
  z[0] = 0.;
  z[1] = 0.;
  CHECK_EQUAL((*f)(z), 0.0);
  z[0] = 10.;
  z[1] = 0.;
  CHECK_EQUAL((*f)(z), 60.0);
  z[0] = 10.;
  z[1] = 5.;
  CHECK_EQUAL((*f)(z), 65.0);
  z[0] = 0.;
  z[1] = 5.;
  CHECK_EQUAL((*f)(z), 5.0);
  // Both coordinates out of bounds
  z[0] = -1.;
  z[1] = -1.;
  CHECK_EQUAL((*f)(z), 0.0);
  z[0] = 11.;
  z[1] = -1.;
  CHECK_EQUAL((*f)(z), 60.0);
  z[0] = 11.;
  z[1] = 6.;
  CHECK_EQUAL((*f)(z), 65.0);
  z[0] = -1.;
  z[1] = 6.;
  CHECK_EQUAL((*f)(z), 5.0);
  // One coordinate out of bounds
  z[0] = 5.5;
  z[1] = -1.;
  CHECK_EQUAL((*f)(z), 33.0);
  z[0] = 5.5;
  z[1] = 6.;
  CHECK_EQUAL((*f)(z), 38.0);
  z[0] = -1.;
  z[1] = 2.5;
  CHECK_EQUAL((*f)(z), 2.5);
  z[0] = 11;
  z[1] = 2.5;
  CHECK_EQUAL((*f)(z), 62.5);
  // Interior points
  z[0] = 5.5;
  z[1] = 2.5;
  CHECK_EQUAL((*f)(z), 35.5);
  z[0] = 2.5;
  z[1] = 2.5;
  CHECK_EQUAL((*f)(z), 17.5);
  z[0] = 7.5;
  z[1] = 2.5;
  CHECK_EQUAL((*f)(z), 47.5);
  z[0] = 5.5;
  z[1] = 1.5;
  CHECK_EQUAL((*f)(z), 34.5);
  z[0] = 5.5;
  z[1] = 3.5;
  CHECK_EQUAL((*f)(z), 36.5);
}

double
f_using_no_params(const Kokkos::View<double*>& x,
                  const Kokkos::View<double*>& p)
{
  return x[1] - x[0];
}

double
f_using_params(const Kokkos::View<double*>& x, const Kokkos::View<double*>& p)
{
  return p[1] * x[1] - p[0] * x[0];
}
