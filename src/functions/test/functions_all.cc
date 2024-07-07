/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//!
#include <iostream>

#include "UnitTest++.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "TestReporterStdout.h"

#include "FunctionConstant.hh"
#include "FunctionSmoothStep.hh"
#include "FunctionTabular.hh"
#include "FunctionPolynomial.hh"
#include "FunctionLinear.hh"
#include "FunctionSeparable.hh"
#include "FunctionStaticHead.hh"
#include "FunctionBilinear.hh"
#include "FunctionComposition.hh"
#include "errors.hh"
#include "HDF5Reader.hh"

using namespace Amanzi;


void
test_values(Function& f, std::vector<double>& values, std::vector<double>& res)
{
  // CPU
  Kokkos::View<double**, Kokkos::HostSpace> x("x", 1,1);
  Kokkos::View<double*, Kokkos::HostSpace> tmp_x("tmp_x", values.size());
  for (int i = 0; i < values.size(); ++i) {
    x(0,0) = values[i];
    CHECK_EQUAL(f(x), res[i]);
    tmp_x(i) = values[i];
  }
  // GPU
  Kokkos::View<double**> x_d("x_d", 1, values.size());
  {
    auto sv = Kokkos::subview(x_d, 0, Kokkos::ALL);
    Kokkos::deep_copy(sv, tmp_x);
  }
  Kokkos::View<double*> res_d("res_d", values.size());
  f.apply(x_d, res_d);
  Kokkos::deep_copy(tmp_x, res_d);
  for (int i = 0; i < values.size(); ++i) CHECK_EQUAL(tmp_x[i], res[i]);
}

void
test_values_multi(Function& f, std::vector<std::vector<double>>& values, std::vector<double>& res)
{
  // CPU
  const int s0 = values.size();
  const int s1 = values[0].size();
  Kokkos::View<double**, Kokkos::HostSpace> x("x", s1, 1);
  Kokkos::DualView<double**> x_d("tmp_x", s1, s0);
  for (int i = 0; i < s0; ++i) {
    for (int j = 0; j < s1; ++j) {
      x(j,0) = values[i][j];
      x_d.view_host()(j, i) = values[i][j];
    }
    CHECK_EQUAL(res[i], f(x));
  }
  // GPU
  Kokkos::deep_copy(x_d.view_device(), x_d.view_host());
  Kokkos::View<double*> res_d("res_d", s0);
  f.apply(x_d.view_device(), res_d);
  Kokkos::View<double*, Kokkos::HostSpace> res_h("res_h", s0);
  Kokkos::deep_copy(res_h, res_d);
  for (int i = 0; i < s0; ++i) CHECK_EQUAL(res[i], res_h[i]);
}


TEST(constant_test)
{
  FunctionConstant f(1.0);
  std::vector<double> val = { 1.0, 3.0 };
  std::vector<double> res = { 1.0, 1.0 };
  test_values(f, val, res);
  auto g = f.Clone();
  test_values(*g, val, res);
}

TEST(smooth_step_test)
{
  FunctionSmoothStep f(1.0, 1.0, 3.0, 0.0);
  std::vector<double> val = { 0.0, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0 };
  std::vector<double> res = { 1.0, 1.0, 0.84375, 0.5, 0.15625, 0.0, 0.0 };
  test_values(f, val, res);
  auto g = f.Clone();
  test_values(*g, val, res);
  CHECK_THROW(FunctionSmoothStep f(3.0, 1.0, 3.0, 0.0), Errors::Message);
}


TEST(tabular_test)
{
  Kokkos::View<double*, Kokkos::HostSpace> x("x", 5);
  x(0) = 0.0;
  x(1) = 1.0;
  x(2) = 3.0;
  x(3) = 3.5;
  x(4) = 3.5;
  Kokkos::View<double*, Kokkos::HostSpace> y("y", 5);
  y(0) = 1.0;
  y(1) = 3.0;
  y(2) = 2.0;
  y(3) = 0.0;
  y(4) = 0.0;
  Kokkos::View<double*, Kokkos::HostSpace> xvec = Kokkos::subview(x, Kokkos::make_pair(0, 4));
  Kokkos::View<double*, Kokkos::HostSpace> yvec = Kokkos::subview(y, Kokkos::make_pair(0, 4));
  int xi = 0;
  FunctionTabular f(xvec, yvec, xi);
  std::vector<double> val = { -1, 0.0, 0.5, 1.0, 2.0, 3.0, 3.25, 3.5, 4.0 };
  std::vector<double> res = { 1.0, 1.0, 2.0, 3.0, 2.5, 2.0, 1.0, 0.0, 0.0 };
  test_values(f, val, res);
  auto g = f.Clone();
  test_values(*g, val, res);

  // Now try with optional form argument
  Kokkos::View<Form_kind*, Kokkos::HostSpace> form("form", 3);
  form(0) = Form_kind::CONSTANT;
  form(1) = Form_kind::LINEAR;
  form(2) = Form_kind::CONSTANT;
  Kokkos::View<Form_kind*, Kokkos::HostSpace> formvec =
    Kokkos::subview(form, Kokkos::make_pair(0, 3));

  FunctionTabular f2(xvec, yvec, xi, formvec);

  std::vector<double> val_2 = { -1, 0.0, 0.5, 1.0, 2.0, 3.0, 3.25, 3.5, 4.0 };
  std::vector<double> res_2 = { 1.0, 1.0, 1.0, 1.0, 2.5, 2.0, 2.0, 2.0, 0.0 };
  test_values(f2, val_2, res_2);
  auto g2 = f2.Clone();
  test_values(*g2, val_2, res_2);

  // Verify the constructor fails with only a single table entry.
  xvec = Kokkos::subview(x, Kokkos::make_pair(0, 1));
  yvec = Kokkos::subview(y, Kokkos::make_pair(0, 1));
  // f = new FunctionTabular(xvec, yvec);
  CHECK_THROW(auto f3 = new FunctionTabular(xvec, yvec, xi), Errors::Message);
  // Verify the constructor fails when the x values are not strictly increasing.
  xvec = Kokkos::subview(x, Kokkos::make_pair(0, 5));
  yvec = Kokkos::subview(y, Kokkos::make_pair(0, 5));
  CHECK_THROW(auto f4 = new FunctionTabular(xvec, yvec, xi), Errors::Message);
  // Verify the constructor fails with different number of x and y values.
  xvec = Kokkos::subview(x, Kokkos::make_pair(0, 4));
  CHECK_THROW(auto f5 = new FunctionTabular(xvec, yvec, xi), Errors::Message);
  // Verify the constructor fails with the wrong number of form values.
  yvec = Kokkos::subview(y, Kokkos::make_pair(0, 4)); // same number of x and y values now
  formvec = Kokkos::subview(formvec, Kokkos::make_pair(0, 2));
  CHECK_THROW(auto f6 = new FunctionTabular(xvec, yvec, xi, formvec), Errors::Message);
}



SUITE(polynomial_test)
{
  TEST(poly1)
  {
    // polynomial with all positive powers
    Kokkos::View<double*, Kokkos::HostSpace> c("c", 2);
    c(0) = -1.0;
    c(1) = 2.0;
    Kokkos::View<int*, Kokkos::HostSpace> p("p", 2);
    p(0) = 3;
    p(1) = 1;
    {
      FunctionPolynomial f(c, p);
      std::vector<double> val = { 2.0, -2.0 };
      std::vector<double> res = { -4.0, 4.0 };
      test_values(f, val, res);

      // same polynomial shifted 2 to the right
      FunctionPolynomial f2(c, p, 2.0);
      std::vector<double> val_2 = { 4.0, 0.0 };
      std::vector<double> res_2 = { -4.0, 4.0 };
      test_values(f2, val_2, res_2);
    }
  }

  TEST(poly2)
  {
    // polynomial with positive powers, constant term, and repeated power
    Kokkos::View<double*, Kokkos::HostSpace> c("c", 4);
    c(0) = -1.0;
    c(1) = 4.0;
    c(2) = 3.0;
    c(3) = -1.0;
    Kokkos::View<int*, Kokkos::HostSpace> p("p", 4);
    p(0) = 3;
    p(1) = 0;
    p(2) = 1;
    p(3) = 1;
    Kokkos::View<double*, Kokkos::HostSpace> cvec = Kokkos::subview(c, Kokkos::make_pair(0, 4));
    Kokkos::View<int*, Kokkos::HostSpace> pvec = Kokkos::subview(p, Kokkos::make_pair(0, 4));
    FunctionPolynomial f(cvec, pvec);
    std::vector<double> val = { 2.0, 0.0, -2.0 };
    std::vector<double> res = { 0.0, 4.0, 8.0 };
    test_values(f, val, res);
  }

  TEST(poly3)
  {
    // polynomial with negative powers only
    // double c[2] = {2.0, -1.0};
    Kokkos::View<double*, Kokkos::HostSpace> c("c", 2);
    c(0) = 2.0;
    c(1) = -1.0;
    // int p[2] = {-3, -2};
    Kokkos::View<int*, Kokkos::HostSpace> p("p", 2);
    p(0) = -3;
    p(1) = -2;
    // std::vector<double> cvec(c, c+2);
    Kokkos::View<double*, Kokkos::HostSpace> cvec = Kokkos::subview(c, Kokkos::make_pair(0, 2));
    // std::vector<int> pvec(p, p+2);
    Kokkos::View<int*, Kokkos::HostSpace> pvec = Kokkos::subview(p, Kokkos::make_pair(0, 2));
    FunctionPolynomial f(cvec, pvec, 1.0);
    // std::vector<double> x(1,3.0);
    std::vector<double> val = { 3.0, 2.0, -1.0 };
    std::vector<double> res = { 0.0, 1.0, -0.5 };
    test_values(f, val, res);
  }

  TEST(poly4)
  {
    // polynomial with positive and negative powers
    Kokkos::View<double*, Kokkos::HostSpace> c("c", 4);
    c(0) = 2.0;
    c(1) = -1.0;
    c(2) = -2.0;
    c(3) = 4.0;
    Kokkos::View<int*, Kokkos::HostSpace> p("p", 4);
    p(0) = 1;
    p(1) = 3;
    p(2) = 0;
    p(3) = -2;
    Kokkos::View<double*, Kokkos::HostSpace> cvec = Kokkos::subview(c, Kokkos::make_pair(0, 4));
    Kokkos::View<int*, Kokkos::HostSpace> pvec = Kokkos::subview(p, Kokkos::make_pair(0, 4));
    FunctionPolynomial f(cvec, pvec);
    std::vector<double> val = { 2.0, 1.0, -1.0, -2.0 };
    std::vector<double> res = { -5.0, 3.0, 1.0, 3.0 };
    test_values(f, val, res);
  }

  TEST(poly_clone)
  {
    // polynomial with positive and negative powers
    Kokkos::View<double*, Kokkos::HostSpace> c("c", 4);
    c(0) = 2.0;
    c(1) = -1.0;
    c(2) = -2.0;
    c(3) = 4.0;
    Kokkos::View<int*, Kokkos::HostSpace> p("p", 4);
    p(0) = 1;
    p(1) = 3;
    p(2) = 0;
    p(3) = -2;
    Kokkos::View<double*, Kokkos::HostSpace> cvec = Kokkos::subview(c, Kokkos::make_pair(0, 4));
    Kokkos::View<int*, Kokkos::HostSpace> pvec = Kokkos::subview(p, Kokkos::make_pair(0, 4));
    FunctionPolynomial f(cvec, pvec);
    std::vector<double> val = { 2.0 };
    std::vector<double> res = { -5.0 };
    test_values(f, val, res);
    auto g = f.Clone();
    test_values(*g, val, res);
  }
}


SUITE(linear_test)
{
  TEST(linear1)
  {
    double y0 = 1.0;
    Kokkos::View<double*, Kokkos::HostSpace> grad("grad", 3);
    grad(0) = 1.0;
    grad(1) = 2.0;
    grad(2) = -3.0;

    // Three variable linear function.
    FunctionLinear f(y0, grad);

    std::vector<std::vector<double>> val = { { 3., 2., 1. }, { 1., -1., 1. } };
    std::vector<double> res = { 5.0, -3.0 };
    test_values_multi(f, val, res);
    auto g = f.Clone();
    test_values_multi(*g, val, res);
  }

  TEST(linear2)
  {
    double y0 = 1.0;
    Kokkos::View<double*, Kokkos::HostSpace> grad("grad", 2);
    grad(0) = 1.0;
    grad(1) = 2.0;

    // Two variable linear function with x0.
    Kokkos::View<double*, Kokkos::HostSpace> x0("x0", 2);
    x0(0) = 1.0;
    x0(1) = 2.0;

    FunctionLinear f2(y0, grad, x0);

    std::vector<std::vector<double>> val = { { 3., 2.}, { 1., -1.} };
    std::vector<double> res = { 3., -5. };
    test_values_multi(f2, val, res);

    // -- check error in case of point smaller than gradient
    Kokkos::View<double**, Kokkos::HostSpace> c("c", 1, 1);
    c(0, 0) = 1.;
    CHECK_THROW(f2(c), Errors::Message);
  }


  TEST(linear3)
  {
    double y0 = 1.0;
    Kokkos::View<double*, Kokkos::HostSpace> grad("grad", 1);
    grad(0) = 1.0;

    FunctionLinear f3(y0, grad);

    std::vector<std::vector<double>> val = { { 3., }, { 1., } };
    std::vector<double> res = { 4., 2. };
    test_values_multi(f3, val, res);
  }

  // Check exception handling.
  TEST(linear_errors)
  {
    double y0 = 1.0;
    Kokkos::View<double*, Kokkos::HostSpace> grad("grad", 0);
    Kokkos::View<double*, Kokkos::HostSpace> x0("x0", 2);
    CHECK_THROW(auto f4 = new FunctionLinear(y0, grad), Errors::Message);
    CHECK_THROW(auto f5 = new FunctionLinear(y0, grad, x0), Errors::Message);

    Kokkos::View<double*, Kokkos::HostSpace> grad2("grad", 1);
    CHECK_THROW(auto f6 = new FunctionLinear(y0, grad, x0), Errors::Message);
  }
}


SUITE(separable_test)
{
  TEST(separable1)
  {
    // First function a smooth step function
    std::unique_ptr<Function> f1 = std::make_unique<FunctionSmoothStep>(1.0, 1.0, 3.0, 0.0);
    // Second function a linear function
    Kokkos::View<double*, Kokkos::HostSpace> grad("grad", 3);
    grad(0) = 1.0;
    grad(1) = 2.0;
    grad(2) = 3.0;
    Kokkos::View<double*, Kokkos::HostSpace> gradvec =
      Kokkos::subview(grad, Kokkos::make_pair(0, 3));
    std::unique_ptr<Function> f2 = std::make_unique<FunctionLinear>(1.0, gradvec);
    FunctionSeparable f(std::move(f1), std::move(f2));
    std::vector<std::vector<double>> val = { { 0., 1., 2., -1. }, { 5.0, 1., 2., -1. } };
    std::vector<double> res = { 3.0, 0.0 };
    test_values_multi(f, val, res);
  }
  TEST(separable2)
  {
    // separable built from separable
    auto fx = std::make_unique<FunctionSmoothStep>(0.0, 1.0, 1.0, 2.0);
    auto fy = std::make_unique<FunctionSmoothStep>(0.0, 1.0, 1.0, 3.0);
    auto fz = std::make_unique<FunctionSmoothStep>(0.0, 1.0, 1.0, 4.0);
    auto fyz = std::make_unique<FunctionSeparable>(std::move(fy), std::move(fz));
    auto fxyz = std::make_unique<FunctionSeparable>(std::move(fx), std::move(fyz));
    std::vector<std::vector<double>> val = {
      { 0., 0., 0. }, { 1., 0., 0. }, { 0., 1., 0. }, { 0., 0., 1. }, { .5, .5, .5 }
    };
    std::vector<double> res = { 1., 2., 3., 4., 7.5 };
    test_values_multi(*fxyz, val, res);
  }
}


TEST(static_head_test)
{
  // txy-linear height function
  double y0 = 1.0;
  Kokkos::View<double*, Kokkos::HostSpace> grad("grad", 4);
  grad(0) = 0.0;
  grad(1) = 2.0;
  grad(2) = 3.0;
  grad(3) = 0.;
  auto h = std::make_unique<FunctionLinear>(y0, grad);
  double patm = -1.0, rho = 0.5, g = 4.0;
  auto f = std::make_unique<FunctionStaticHead>(patm, rho, g, std::move(h), 3);
  std::vector<std::vector<double>> val = { { 0., 1., .5, 2. }, { 1., 2., 1., 1. } };
  std::vector<double> res = {
    patm + rho * g * (y0 + grad[1] * val[0][1] + grad[2] * val[0][2] - val[0][3]),
    (patm + rho * g * (y0 + grad[1] * val[1][1] + grad[2] * val[1][2] - val[1][3]))
  };
  auto ff = f->Clone();
  test_values_multi(*ff, val, res);
}


TEST(bilinear_test)
{
  std::string filename = "test/bilinear.h5";
  HDF5Reader reader(filename);
  std::string row_name = "times";
  // std::vector<double> vec_x;
  // HostSpace test
  Kokkos::View<double*, Kokkos::HostSpace> vec_x;
  int xi = 0;
  reader.ReadData(row_name, vec_x);
  std::string col_name = "x";
  // std::vector<double> vec_y;
  Kokkos::View<double*, Kokkos::HostSpace> vec_y;
  int yi = 1;
  reader.ReadData(col_name, vec_y);
  std::string v_name = "values";

  Kokkos::View<double**, Kokkos::HostSpace> mat_v;
  reader.ReadMatData(v_name, mat_v);
  auto f = std::make_unique<FunctionBilinear>(vec_x, vec_y, mat_v, xi, yi);
  // Corners
  // std::vector<double> z(2,0.);
  std::vector<std::vector<double>> val = { { 0., 0. },   { 10., 0. },  { 10., 5. },  { 0., 5. },
                                           { -1., -1 },  { 11., -1. }, { 11., 6. },  { -1., 6. },
                                           { 5.5, -1. }, { 5.5, 6. },  { -1., 2.5 }, { 11, 2.5 },
                                           { 5.5, 2.5 }, { 2.5, 2.5 }, { 7.5, 2.5 }, { 5.5, 1.5 },
                                           { 5.5, 3.5 } };
  std::vector<double> res = { 0.0,  60.0, 65.0, 5.0,  0.0,  60.0, 65.0, 5.0, 33.,
                              38.0, 2.5,  62.5, 35.5, 17.5, 47.5, 34.5, 36.5 };
  test_values_multi(*f, val, res);
}


TEST(composition_test)
{
  double y0 = 1.0;
  Kokkos::View<double*, Kokkos::HostSpace> grad("grad", 3);
  grad(0) = 1.0;
  grad(1) = 2.0;
  grad(2) = -3.0;
  Kokkos::View<double*, Kokkos::HostSpace> gradvec = Kokkos::subview(grad, Kokkos::make_pair(0, 3));
  // Three variable linear function.
  auto f2 = std::make_unique<FunctionLinear>(y0, gradvec);


  Kokkos::View<double*, Kokkos::HostSpace> x("x", 3);
  x(0) = -3.0;
  x(1) = 0.0;
  x(2) = 5.0;
  Kokkos::View<double*, Kokkos::HostSpace> y("y", 3);
  y(0) = 1.0;
  y(1) = 3.0;
  y(2) = 2.0;
  int xi = 0;
  auto f1 = std::make_unique<FunctionTabular>(x, y, xi);

  FunctionComposition f(std::move(f1), std::move(f2));

  std::vector<std::vector<double>> val = { { 3., 2., 1. }, { 1., -1., 1. } };
  std::vector<double> res = { 2., 1. };
  test_values_multi(f, val, res);
}

