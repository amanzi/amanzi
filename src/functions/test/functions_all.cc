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
#include "errors.hh"
#include "HDF5Reader.hh"

using namespace Amanzi;


void test_values(Function* f, std::vector<double>& values,std::vector<double>& res){ 
  // CPU
  Kokkos::View<double*,Kokkos::HostSpace> x("x", 1);
  Kokkos::View<double*,Kokkos::HostSpace> tmp_x("tmp_x",values.size()); 
  for(int i = 0 ; i < values.size() ; ++i){
    x(0) = values[i]; 
    CHECK_EQUAL((*f)(x), res[i]); 
    tmp_x(i) = values[i];  
  }
  // GPU
  Kokkos::View<double**> x_d("x_d",1,values.size());  
  auto sv = Kokkos::subview(x_d,0,Kokkos::ALL);
  Kokkos::deep_copy(sv,tmp_x); 
  Kokkos::View<double*> res_d("res_d",values.size());
  f->apply(x_d,res_d);
  Kokkos::deep_copy(tmp_x,res_d);
  for(int i = 0 ; i < values.size(); ++i)
    CHECK_EQUAL(tmp_x[i],res[i]);   
}

void test_values_multi(Function* f, std::vector<std::vector<double>>& values, std::vector<double> & res){ 
  // CPU
  const int s0 = values.size(); 
  const int s1 = values[0].size(); 
  Kokkos::View<double*,Kokkos::HostSpace> x("x", s1);
  Kokkos::DualView<double**> x_d("tmp_x",s1,s0); 
  for(int i = 0 ; i < s0 ; ++i){
    for(int j = 0 ; j < s1 ; ++j){
      x(j) = values[i][j]; 
      x_d.view_host()(j,i) = values[i][j]; 
    }
    CHECK_EQUAL((*f)(x), res[i]); 
  }
  // GPU
  Kokkos::deep_copy(x_d.view_device(),x_d.view_host());  
  Kokkos::View<double*> res_d("res_d",s0);
  f->apply(x_d.view_device(),res_d);
  Kokkos::View<double*,Kokkos::HostSpace> res_h("res_h",s0); 
  Kokkos::deep_copy(res_h,res_d);
  for(int i = 0 ; i < s0; ++i)
    CHECK_EQUAL(res_h[i],res[i]);   
}

int
main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize();
  int status = UnitTest::RunAllTests();
  Kokkos::finalize();
  return status;
}

TEST(constant_test)
{
  Function* f = new FunctionConstant(1.0);
  std::vector<double> val = {1.0,3.0};
  std::vector<double> res = {1.0,1.0};  
  test_values(f,val,res); 
  Function* g = f->Clone();
  delete f;
  test_values(g,val,res); 
  delete g;
}

TEST(smooth_step_test)
{
  Function* f = new FunctionSmoothStep(1.0, 1.0, 3.0, 0.0);
  std::vector<double> val = {0.0,1.0,1.5,2.0,2.5,3.0,4.0}; 
  std::vector<double> res = {1.0,1.0,0.84375,0.5,0.15625,0.0,0.0}; 
  test_values(f,val,res); 
  Function* g = f->Clone();
  delete f;
  test_values(g,val,res); 
  delete g;
  CHECK_THROW(FunctionSmoothStep f(3.0, 1.0, 3.0, 0.0), Errors::Message);
}


TEST(tabular_test)
{
  Kokkos::View<double*,Kokkos::HostSpace> x("x", 5);
  x(0) = 0.0;
  x(1) = 1.0;
  x(2) = 3.0;
  x(3) = 3.5;
  x(4) = 3.5;
  Kokkos::View<double*,Kokkos::HostSpace> y("y", 5);
  y(0) = 1.0;
  y(1) = 3.0;
  y(2) = 2.0;
  y(3) = 0.0;
  y(4) = 0.0;
  Kokkos::View<double*,Kokkos::HostSpace> xvec = Kokkos::subview(x, Kokkos::make_pair(0, 4));
  Kokkos::View<double*,Kokkos::HostSpace> yvec = Kokkos::subview(y, Kokkos::make_pair(0, 4));
  int xi = 0;
  Function* f = new FunctionTabular(xvec, yvec, xi);
  std::vector<double> val = {-1,0.0,0.5,1.0,2.0,3.0,3.25,3.5,4.0}; 
  std::vector<double> res = {1.0,1.0,2.0,3.0,2.5,2.0,1.0,0.0,0.0}; 
  test_values(f,val,res); 
  Function* g = f->Clone();
  delete f;
  test_values(g,val,res); 
  delete g;

  // Now try with optional form argument
  Kokkos::View<FunctionTabular::Form*,Kokkos::HostSpace> form("form", 3);
  form(0) = FunctionTabular::CONSTANT;
  form(1) = FunctionTabular::LINEAR;
  form(2) = FunctionTabular::CONSTANT;
  Kokkos::View<FunctionTabular::Form*,Kokkos::HostSpace> formvec =
    Kokkos::subview(form, Kokkos::make_pair(0, 3));

  f = new FunctionTabular(xvec, yvec, xi, formvec);

  std::vector<double> val_2 = {-1,0.0,0.5,1.0,2.0,3.0,3.25,3.5,4.0}; 
  std::vector<double> res_2 = {1.0,1.0,1.0,1.0,2.5,2.0,2.0,2.0,0.0}; 
  test_values(f,val_2,res_2); 
  g = f->Clone();
  delete f;
  test_values(g,val_2,res_2); 
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
    Kokkos::View<double*,Kokkos::HostSpace> c("c", 2);
    c(0) = -1.0;
    c(1) = 2.0;
    Kokkos::View<int*,Kokkos::HostSpace> p("p", 2);
    p(0) = 3;
    p(1) = 1;
    Kokkos::View<double*,Kokkos::HostSpace> cvec = Kokkos::subview(c, Kokkos::make_pair(0, 2));
    Kokkos::View<int*,Kokkos::HostSpace> pvec = Kokkos::subview(p, Kokkos::make_pair(0, 2));
    Function* f = new FunctionPolynomial(cvec, pvec);

    std::vector<double> val = {2.0,-2.0}; 
    std::vector<double> res = {-4.0,4.0}; 
    test_values(f,val,res); 
    delete f;
    // same polynomial shifted 2 to the right
    f = new FunctionPolynomial(cvec, pvec, 2.0);
    std::vector<double> val_2 = {4.0,0.0}; 
    std::vector<double> res_2 = {-4.0,4.0}; 
    test_values(f,val_2,res_2); 
    delete f;
  }
  TEST(poly2)
  {
    // polynomial with positive powers, constant term, and repeated power
    Kokkos::View<double*,Kokkos::HostSpace> c("c", 4);
    c(0) = -1.0;
    c(1) = 4.0;
    c(2) = 3.0;
    c(3) = -1.0;
    Kokkos::View<int*,Kokkos::HostSpace> p("p", 4);
    p(0) = 3;
    p(1) = 0;
    p(2) = 1;
    p(3) = 1;
    Kokkos::View<double*,Kokkos::HostSpace> cvec = Kokkos::subview(c, Kokkos::make_pair(0, 4));
    Kokkos::View<int*,Kokkos::HostSpace> pvec = Kokkos::subview(p, Kokkos::make_pair(0, 4));
    Function* f = new FunctionPolynomial(cvec, pvec);
    std::vector<double> val = {2.0,0.0,-2.0}; 
    std::vector<double> res = {0.0,4.0,8.0};
    test_values(f,val,res);  
    delete f;
  }
  TEST(poly3)
  {
    // polynomial with negative powers only
    // double c[2] = {2.0, -1.0};
    Kokkos::View<double*,Kokkos::HostSpace> c("c", 2);
    c(0) = 2.0;
    c(1) = -1.0;
    // int p[2] = {-3, -2};
    Kokkos::View<int*,Kokkos::HostSpace> p("p", 2);
    p(0) = -3;
    p(1) = -2;
    // std::vector<double> cvec(c, c+2);
    Kokkos::View<double*,Kokkos::HostSpace> cvec = Kokkos::subview(c, Kokkos::make_pair(0, 2));
    // std::vector<int> pvec(p, p+2);
    Kokkos::View<int*,Kokkos::HostSpace> pvec = Kokkos::subview(p, Kokkos::make_pair(0, 2));
    Function* f = new FunctionPolynomial(cvec, pvec, 1.0);
    // std::vector<double> x(1,3.0);   
    std::vector<double> val = {3.0,2.0,-1.0}; 
    std::vector<double> res = {0.0,1.0,-0.5};
    test_values(f,val,res);  
    delete f;
  }
  TEST(poly4)
  {
    // polynomial with positive and negative powers
    Kokkos::View<double*,Kokkos::HostSpace> c("c", 4);
    c(0) = 2.0;
    c(1) = -1.0;
    c(2) = -2.0;
    c(3) = 4.0;
    Kokkos::View<int*,Kokkos::HostSpace> p("p", 4);
    p(0) = 1;
    p(1) = 3;
    p(2) = 0;
    p(3) = -2;
    Kokkos::View<double*,Kokkos::HostSpace> cvec = Kokkos::subview(c, Kokkos::make_pair(0, 4));
    Kokkos::View<int*,Kokkos::HostSpace> pvec = Kokkos::subview(p, Kokkos::make_pair(0, 4));
    Function* f = new FunctionPolynomial(cvec, pvec);
    std::vector<double> val = {2.0,1.0,-1.0,-2.0};
    std::vector<double> res = {-5.0,3.0,1.0,3.0};  
    test_values(f,val,res); 
    delete f;
  }
  TEST(poly_clone)
  {
    // polynomial with positive and negative powers
    Kokkos::View<double*,Kokkos::HostSpace> c("c", 4);
    c(0) = 2.0;
    c(1) = -1.0;
    c(2) = -2.0;
    c(3) = 4.0;
    Kokkos::View<int*,Kokkos::HostSpace> p("p", 4);
    p(0) = 1;
    p(1) = 3;
    p(2) = 0;
    p(3) = -2;
    Kokkos::View<double*,Kokkos::HostSpace> cvec = Kokkos::subview(c, Kokkos::make_pair(0, 4));
    Kokkos::View<int*,Kokkos::HostSpace> pvec = Kokkos::subview(p, Kokkos::make_pair(0, 4));
    FunctionPolynomial* f = new FunctionPolynomial(cvec, pvec);
    std::vector<double> val = {2.0};  
    std::vector<double> res = {-5.0}; 
    test_values(f,val,res); 
    FunctionPolynomial* g = f->Clone();
    delete f;
    test_values(g,val,res); 
    delete g;
  }
}


TEST(linear_test)
{
  double y0 = 1.0;
  Kokkos::View<double*,Kokkos::HostSpace> grad("grad", 3);
  grad(0) = 1.0;
  grad(1) = 2.0;
  grad(2) = -3.0;
  Kokkos::View<double*,Kokkos::HostSpace> gradvec =
    Kokkos::subview(grad, Kokkos::make_pair(0, 3));
  Function* f;
  // Three variable linear function.
  f = new FunctionLinear(y0, gradvec);

  std::vector<std::vector<double>> val = {{3.,2.,1.},{1.,-1.,1.}}; 
  std::vector<double> res = {5.0,-3.0}; 
  test_values_multi(f,val,res); 
  Function* g = f->Clone();
  delete f;
  test_values_multi(g,val,res);
  delete g;
  // Two variable linear function with x0.
  Kokkos::View<double*,Kokkos::HostSpace> x0("x0", 3);
  x0(0) = 1.0;
  x0(1) = 2.0;
  x0(2) = 3.0;
  Kokkos::View<double*,Kokkos::HostSpace> x0vec = Kokkos::subview(x0, Kokkos::make_pair(0, 2));
  gradvec = Kokkos::subview(grad, Kokkos::make_pair(0, 2));
  f = new FunctionLinear(y0, gradvec, x0vec);
  std::vector<double> res_2 = {3.,-5.}; 
  test_values_multi(f,val,res_2); 
  // -- check error in case of point smaller than gradient
  Kokkos::View<double*,Kokkos::HostSpace> c("c", 1);
  c(0) = 1.;
  CHECK_THROW((*f)(c), Errors::Message);
  delete f;

  // Single variable linear function.
  gradvec = Kokkos::subview(grad, Kokkos::make_pair(0, 1));
  f = new FunctionLinear(y0, gradvec);
  std::vector<double> res_3 = {4.,2.}; 
  test_values_multi(f,val,res_3); 
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
    Kokkos::View<double*,Kokkos::HostSpace> grad("grad", 3);
    grad(0) = 1.0;
    grad(1) = 2.0;
    grad(2) = 3.0;
    Kokkos::View<double*,Kokkos::HostSpace> gradvec =
      Kokkos::subview(grad, Kokkos::make_pair(0, 3));
    Function* f2 = new FunctionLinear(1.0, gradvec);
    Function* f = new FunctionSeparable(f1, *f2);
    delete f2;
    std::vector<std::vector<double>> val = {{0.,1.,2.,-1.},{5.0,1.,2.,-1.}};
    std::vector<double> res = {3.0,0.0};
    test_values_multi(f,val,res);  
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
    std::vector<std::vector<double>> val = {{0.,0.,0.},{1.,0.,0.},
      {0.,1.,0.},{0.,0.,1.},{.5,.5,.5}}; 
    std::vector<double> res = {1.,2.,3.,4.,7.5}; 
    test_values_multi(fxyz,val,res); 
    delete fxyz;
  }
}

TEST(static_head_test)
{
  // txy-linear height function
  double y0 = 1.0;
  Kokkos::View<double*,Kokkos::HostSpace> grad("grad", 3);
  grad(0) = 0.0;
  grad(1) = 2.0;
  grad(2) = 3.0;
  Kokkos::View<double*,Kokkos::HostSpace> gradvec =
    Kokkos::subview(grad, Kokkos::make_pair(0, 3));
  Function* h = new FunctionLinear(y0, gradvec);
  double patm = -1.0, rho = 0.5, g = 4.0;
  Function* f = new FunctionStaticHead(patm, rho, g, *h, 3);
  std::vector<std::vector<double>> val = {{0.,1.,.5,2.},{1.,2.,1.,1.}};
  std::vector<double> res = {patm + rho * g * (y0 + grad[1] * val[0][1] + grad[2] * val[0][2] - val[0][3]),(patm + rho * g * (y0 + grad[1] * val[1][1] + grad[2] * val[1][2] - val[1][3]))}; 
  Function* ff = f->Clone();
  delete f;
  test_values_multi(ff,val,res); 
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
  // HostSpace test 
  Kokkos::View<double*,Kokkos::HostSpace> vec_x;
  int xi = 0;
  reader.ReadData(row_name, vec_x);
  std::string col_name = "x";
  // std::vector<double> vec_y;
  Kokkos::View<double*,Kokkos::HostSpace> vec_y;
  int yi = 1;
  reader.ReadData(col_name, vec_y);
  std::string v_name = "values";

  Kokkos::View<double**,Kokkos::HostSpace> mat_v;
  reader.ReadMatData(v_name, mat_v);
  Function* f = new FunctionBilinear(vec_x, vec_y, mat_v, xi, yi);
  // Corners
  // std::vector<double> z(2,0.);
  std::vector<std::vector<double>> val = 
    {{0.,0.},{10.,0.},{10.,5.},{0.,5.},
     {-1.,-1},{11.,-1.},{11.,6.},{-1.,6.},
     {5.5,-1.},{5.5,6.},{-1.,2.5},{11,2.5},
     {5.5,2.5},{2.5,2.5},{7.5,2.5},{5.5,1.5},{5.5,3.5}}; 
  std::vector<double> res = 
    {0.0,60.0,65.0,5.0,
     0.0,60.0,65.0,5.0,
     33.,38.0,2.5,62.5,
     35.5,17.5,47.5,34.5,36.5}; 
  test_values_multi(f,val,res); 
}

double
f_using_no_params(const Kokkos::View<double*,Kokkos::HostSpace>& x,
                  const Kokkos::View<double*,Kokkos::HostSpace>& p)
{
  return x[1] - x[0];
}

double
f_using_params(const Kokkos::View<double*,Kokkos::HostSpace>& x, const Kokkos::View<double*,Kokkos::HostSpace>& p)
{
  return p[1] * x[1] - p[0] * x[0];
}
