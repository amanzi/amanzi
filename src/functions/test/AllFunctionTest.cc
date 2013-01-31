#include "UnitTest++.h"
#include "TestReporterStdout.h"

#include "ConstantFunction.hh"
#include "SmoothStepFunction.hh"
#include "TabularFunction.hh"
#include "PolynomialFunction.hh"
#include "LinearFunction.hh"
#include "SeparableFunction.hh"
#include "PointerFunction.hh"
#include "StaticHeadFunction.hh"
#include "errors.hh"

#include <iostream>

using namespace Amanzi;

int main (int argc, char *argv[])
{
    return UnitTest::RunAllTests ();
}

TEST(constant_test)
{
  Function *f = new ConstantFunction(1.0);
  double x = 3.0;
  CHECK_EQUAL((*f)(&x), 1.0);
  Function *g = f->Clone();
  delete f;
  CHECK_EQUAL((*g)(&x), 1.0);
  delete g;
}

TEST(smooth_step_test)
{
  Function *f = new SmoothStepFunction(1.0, 1.0, 3.0, 0.0);
  double x;
  x = 0.0;
  CHECK_EQUAL((*f)(&x), 1.0);
  x = 1.0;
  CHECK_EQUAL((*f)(&x), 1.0);
  x = 1.5;
  CHECK_EQUAL((*f)(&x), 0.84375);
  x = 2.0;
  CHECK_EQUAL((*f)(&x), 0.5);
  x = 2.5;
  CHECK_EQUAL((*f)(&x), 0.15625);
  x = 3.0;
  CHECK_EQUAL((*f)(&x), 0.0);
  x = 4.0;
  CHECK_EQUAL((*f)(&x), 0.0);
  Function *g = f->Clone();
  delete f;
  x = 2.0;
  CHECK_EQUAL((*g)(&x), 0.5);
  delete g;
  CHECK_THROW(SmoothStepFunction f(3.0, 1.0, 3.0, 0.0), Errors::Message);
}
      
TEST(tabular_test)
{
  double x[5] = {0.0, 1.0, 3.0, 3.5, 3.5};
  double y[5] = {1.0, 3.0, 2.0, 0.0, 0.0};
  std::vector<double> xvec(x, x+4);
  std::vector<double> yvec(y, y+4);
  Function *f = new TabularFunction(xvec, yvec);
  double z;
  z = -1.0;
  CHECK_EQUAL((*f)(&z), 1.0);
  z = 0.0;
  CHECK_EQUAL((*f)(&z), 1.0);
  z = 0.5;
  CHECK_EQUAL((*f)(&z), 2.0);
  z = 1.0;
  CHECK_EQUAL((*f)(&z), 3.0);
  z = 2.0;
  CHECK_EQUAL((*f)(&z), 2.5);
  z = 3.0;
  CHECK_EQUAL((*f)(&z), 2.0);
  z = 3.25;
  CHECK_EQUAL((*f)(&z), 1.0);
  z = 3.5;
  CHECK_EQUAL((*f)(&z), 0.0);
  z = 4.0;
  CHECK_EQUAL((*f)(&z), 0.0);
  Function *g = f->Clone();
  delete f;
  z = 3.0;
  CHECK_EQUAL((*g)(&z), 2.0);
  delete g;
  // Now try with optional form argument
  TabularFunction::Form form[3] = {TabularFunction::CONSTANT,
      TabularFunction::LINEAR, TabularFunction::CONSTANT};
  std::vector<TabularFunction::Form> formvec(form, form+3);
  f = new TabularFunction(xvec, yvec, formvec);
  z = -1.0;
  CHECK_EQUAL((*f)(&z), 1.0);
  z = 0.0;
  CHECK_EQUAL((*f)(&z), 1.0);
  z = 0.5;
  CHECK_EQUAL((*f)(&z), 1.0);
  z = 1.0;
  CHECK_EQUAL((*f)(&z), 1.0);
  z = 2.0;
  CHECK_EQUAL((*f)(&z), 2.5);
  z = 3.0;
  CHECK_EQUAL((*f)(&z), 2.0);
  z = 3.25;
  CHECK_EQUAL((*f)(&z), 2.0);
  z = 3.5;
  CHECK_EQUAL((*f)(&z), 2.0);
  z = 4.0;
  CHECK_EQUAL((*f)(&z), 0.0);
  g = f->Clone();
  delete f;
  z = 3.0;
  CHECK_EQUAL((*g)(&z), 2.0);
  delete g;
  // Verify the constructor fails with only a single table entry.
  xvec.assign(x, x+1);
  yvec.assign(y, y+1);
  //f = new TabularFunction(xvec, yvec);
  CHECK_THROW(f = new TabularFunction(xvec, yvec), Errors::Message);
  // Verify the constructor fails when the x values are not strictly increasing.
  xvec.assign(x, x+5);
  yvec.assign(y, y+5);
  //f = new TabularFunction(xvec, yvec);
  CHECK_THROW(f = new TabularFunction(xvec, yvec), Errors::Message);
  // Verify the constructor fails with different number of x and y values.
  xvec.pop_back();
  //f = new TabularFunction(xvec, yvec);
  CHECK_THROW(f = new TabularFunction(xvec, yvec), Errors::Message);
  // Verify the constructor fails with the wrong number of form values.
  yvec.pop_back(); // same number of x and y values now
  formvec.pop_back();
  //f = new TabularFunction(xvec, yvec, formvec);
  CHECK_THROW(f = new TabularFunction(xvec, yvec, formvec), Errors::Message);
}

SUITE(polynomial_test) {
  TEST(poly1)
  {
    // polynomial with all positive powers
    double c[2] = {-1.0, 2.0};
    int p[2] = {3, 1};
    std::vector<double> cvec(c, c+2);
    std::vector<int> pvec(p, p+2);
    Function *f = new PolynomialFunction(cvec, pvec);
    double x = 2.0;
    CHECK_EQUAL((*f)(&x), -4.0);
    x = -2.0;
    CHECK_EQUAL((*f)(&x), 4.0);
    delete f;
    // same polynomial shifted 2 to the right
    f = new PolynomialFunction(cvec, pvec, 2.0);
    x = 4.0;
    CHECK_EQUAL((*f)(&x), -4.0);
    x = 0.0;
    CHECK_EQUAL((*f)(&x), 4.0);
    delete f;
  }
  TEST(poly2)
  {
    // polynomial with positive powers, constant term, and repeated power
    double c[4] = {-1.0, 4.0, 3.0, -1.0};
    int p[4] = {3, 0, 1, 1};
    std::vector<double> cvec(c, c+4);
    std::vector<int> pvec(p, p+4);
    Function *f = new PolynomialFunction(cvec, pvec);
    double x = 2.0;
    CHECK_EQUAL((*f)(&x), 0.0);
    x = 0.0;
    CHECK_EQUAL((*f)(&x), 4.0);
    x = -2.0;
    CHECK_EQUAL((*f)(&x), 8.0);
    delete f;
  }
  TEST(poly3)
  {
    // polynomial with negative powers only
    double c[2] = {2.0, -1.0};
    int p[2] = {-3, -2};
    std::vector<double> cvec(c, c+2);
    std::vector<int> pvec(p, p+2);
    Function *f = new PolynomialFunction(cvec, pvec, 1.0);
    double x = 3.0;
    CHECK_EQUAL((*f)(&x), 0.0);
    x = 2.0;
    CHECK_EQUAL((*f)(&x), 1.0);
    x = -1.0;
    CHECK_EQUAL((*f)(&x), -0.5);
    delete f;
  }
  TEST(poly4)
  {
    // polynomial with positive and negative powers
    double c[4] = {2.0, -1.0, -2.0, 4.0};
    int p[4] = {1, 3, 0, -2};
    std::vector<double> cvec(c, c+4);
    std::vector<int> pvec(p, p+4);
    Function *f = new PolynomialFunction(cvec, pvec);
    double x = 2.0;
    CHECK_EQUAL((*f)(&x), -5.0);
    x = 1.0;
    CHECK_EQUAL((*f)(&x), 3.0);
    x = -1.0;
    CHECK_EQUAL((*f)(&x), 1.0);
    x = -2.0;
    CHECK_EQUAL((*f)(&x), 3.0);
    delete f;
  }
  TEST(poly_clone)
  {
    // polynomial with positive and negative powers
    double c[4] = {2.0, -1.0, -2.0, 4.0};
    int p[4] = {1, 3, 0, -2};
    std::vector<double> cvec(c, c+4);
    std::vector<int> pvec(p, p+4);
    PolynomialFunction *f = new PolynomialFunction(cvec, pvec);
    double x = 2.0;
    CHECK_EQUAL((*f)(&x), -5.0);
    PolynomialFunction *g = f->Clone();
    delete f;
    CHECK_EQUAL((*g)(&x), -5.0);
    delete g;
  }
}

TEST(linear_test)
{
  double y0 = 1.0;
  double grad[3] = {1.0, 2.0, -3.0};
  std::vector<double> gradvec(grad, grad+3);
  Function *f;
  // Three variable linear function.
  f = new LinearFunction(y0, gradvec);
  double a[3] = {3.0, 2.0, 1.0};
  double b[3] = {1.0, -1.0, 1.0};
  CHECK_EQUAL(5.0, (*f)(a));
  CHECK_EQUAL(-3.0, (*f)(b));
  Function *g = f->Clone();
  delete f;
  CHECK_EQUAL(-3.0, (*g)(b));
  delete g;
  // Two variable linear function with x0.
  double x0[3] = { 1.0, 2.0, 3.0};
  std::vector<double> x0vec(x0, x0+2);
  gradvec.pop_back();
  f = new LinearFunction(y0, gradvec, x0vec);
  CHECK_EQUAL(3.0, (*f)(a));
  CHECK_EQUAL(-5.0, (*f)(b));
  delete f;
  // Single variable linear function.
  gradvec.pop_back();
  f = new LinearFunction(y0, gradvec);
  CHECK_EQUAL(4.0, (*f)(a));
  CHECK_EQUAL(2.0, (*f)(b));
  delete f;
  // Check exception handling.
  gradvec.pop_back(); // 0-sized now
  //f = new LinearFunction(y0, gradvec);
  CHECK_THROW(f = new LinearFunction(y0, gradvec), Errors::Message);
  //f = new LinearFunction(y0, gradvec, x0vec);
  CHECK_THROW(f = new LinearFunction(y0, gradvec, x0vec), Errors::Message);
  gradvec.push_back(1.0); // 1-sized, but x0vec is 2-sized
  //f = new LinearFunction(y0, gradvec, x0vec);
  CHECK_THROW(f = new LinearFunction(y0, gradvec, x0vec), Errors::Message);
}

SUITE(separable_test) {
  TEST(separable1)
  {
    // First function a smooth step function
    SmoothStepFunction f1(1.0, 1.0, 3.0, 0.0);
    // Second function a linear function
    double grad[3] = {1.0, 2.0, 3.0};
    std::vector<double> gradvec(grad, grad+3);
    Function *f2 = new LinearFunction(1.0, gradvec);
    Function *f = new SeparableFunction(f1, *f2);
    delete f2;
    double x[4] = {0.0, 1.0, 2.0, -1.0}; // t, x, y, z
    CHECK_EQUAL(3.0, (*f)(x));
    x[0] = 5.0;
    CHECK_EQUAL(0.0, (*f)(x));
    delete f;
  }

  TEST(separable2)
  {
    // separable built from separable
    Function *fx = new SmoothStepFunction(0.0, 1.0, 1.0, 2.0);
    Function *fy = new SmoothStepFunction(0.0, 1.0, 1.0, 3.0);
    Function *fz = new SmoothStepFunction(0.0, 1.0, 1.0, 4.0);
    Function *fyz = new SeparableFunction(*fy, *fz);
    Function *fxyz = new SeparableFunction(*fx, *fyz);
    delete fx, fy, fz, fyz;
    double x0[3] = {0.0, 0.0, 0.0};
    double x1[3] = {1.0, 0.0, 0.0};
    double x2[3] = {0.0, 1.0, 0.0};
    double x3[3] = {0.0, 0.0, 1.0};
    double x4[3] = {0.5, 0.5, 0.5};
    CHECK_EQUAL(1.0, (*fxyz)(x0));
    CHECK_EQUAL(2.0, (*fxyz)(x1));
    CHECK_EQUAL(3.0, (*fxyz)(x2));
    CHECK_EQUAL(4.0, (*fxyz)(x3));
    CHECK_EQUAL(7.5, (*fxyz)(x4));
    delete fxyz;
  }
}

TEST(static_head_test)
{
  // txy-linear height function
  double y0 = 1.0, grad[3] = {0.0, 2.0, 3.0};
  std::vector<double> gradvec(grad, grad+3);
  Function *h = new LinearFunction(y0, gradvec);
  double patm = -1.0, rho = 0.5, g = 4.0;
  Function *f = new StaticHeadFunction(patm, rho, g, *h, 3);
  double a[4] = {0.0, 1.0, 0.5, 2.0}; // t, x, y, z
  CHECK_EQUAL(patm+rho*g*(y0+grad[1]*a[1]+grad[2]*a[2]-a[3]), (*f)(a));
  double b[4] = {1.0, 2.0, 1.0, 1.0}; // t, x, y, z
  CHECK_EQUAL(patm+rho*g*(y0+grad[1]*b[1]+grad[2]*b[2]-b[3]), (*f)(b));
  Function *ff = f->Clone();
  double val = (*f)(b);
  delete f;
  CHECK_EQUAL(val, (*ff)(b));
}

//SUITE(separable_test) {
//  TEST(separable1)
//  {
//    // First function a smooth step function
//    std::auto_ptr<Function> f1(new SmoothStepFunction(1.0, 1.0, 3.0, 0.0));
//    // Second function a linear function
//    double grad[3] = {1.0, 2.0, 3.0};
//    std::vector<double> gradvec(grad, grad+3);
//    std::auto_ptr<Function> f2(new LinearFunction(1.0, gradvec));
//    Function *f = new SeparableFunction(f1, f2);
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
//    std::auto_ptr<Function> fx(new SmoothStepFunction(0.0, 1.0, 1.0, 2.0));
//    std::auto_ptr<Function> fy(new SmoothStepFunction(0.0, 1.0, 1.0, 3.0));
//    std::auto_ptr<Function> fz(new SmoothStepFunction(0.0, 1.0, 1.0, 4.0));
//    std::auto_ptr<Function> fyz(new SeparableFunction(fy, fz));
//    Function *fxyz = new SeparableFunction(fx, fyz);
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

double f_using_no_params (const double *x, const double *p)
{
  return x[1] - x[0];
}

double f_using_params (const double *x, const double *p)
{
  return p[1]*x[1] - p[0]*x[0];
}

TEST(pointer_test)
{
  Function *f = new PointerFunction(&f_using_no_params);
  double x[2] = {2.0, 1.0};
  CHECK_EQUAL(f_using_no_params(x,0), (*f)(x));
  delete f;
  double p[2] = {-1.0, 2.0};
  std::vector<double> pvec(p,p+2);
  f = new PointerFunction(&f_using_params, pvec);
  CHECK_EQUAL(f_using_params(x,p), (*f)(x));
  Function *g = f->Clone();
  delete f;
  CHECK_EQUAL(f_using_params(x,p), (*g)(x));
  delete g;
}
