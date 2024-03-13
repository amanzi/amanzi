/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Operators

  Test for a nonlinear coupled diffusion problem.

  - div( k00(u,v) grad(u) ) = f0
  - div( k11(u,v) grad(v) ) = f1

  u = a x^2 + y^2
  v = sin^2( pi/4 * (x+y) )

  k0 = exp(uv)
  k1 = exp((u+v)/2)

  Note that k0 and k1 are both positive, monotonically increasing functions of
  both u and v.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_NONLINEAR_00_HH_
#define AMANZI_OPERATOR_ANALYTIC_NONLINEAR_00_HH_

#include "AnalyticNonlinearCoupledBase.hh"

class AnalyticNonlinearCoupled00 : public AnalyticNonlinearCoupledBase {
 public:
  AnalyticNonlinearCoupled00(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh)
    : AnalyticNonlinearCoupledBase(mesh){};

  // analytic solution for diffusion problem with gravity
  bool isBlock(int i, int j) { return i == j ? true : false; }

  double ScalarCoefficient00(double u, double v) { return exp(u * v); }

  double ScalarCoefficient11(double u, double v) { return exp((u + v) / 2.); }

  double DScalarCoefficient00D0(double u, double v)
  {
    // if (exp(u*v) * v < 0) {
    //   std::cout << "negative deriv" << std::endl;
    // }
    return exp(u * v) * v;
  }
  double DScalarCoefficient00D1(double u, double v)
  {
    // if (exp(u*v) * u < 0) {
    //   std::cout << "negative deriv" << std::endl;
    // }
    return exp(u * v) * u;
  }

  double DScalarCoefficient11D0(double u, double v)
  {
    // if (exp((u+v)/2.) * 0.5 < 0) {
    //   std::cout << "negative deriv" << std::endl;
    // }
    return exp((u + v) / 2.) * 0.5;
  }
  double DScalarCoefficient11D1(double u, double v)
  {
    // if (exp((u+v)/2.) * 0.5 < 0) {
    //   std::cout << "negative deriv" << std::endl;
    // }
    return exp((u + v) / 2.) * 0.5;
  }

  double exact0(const Amanzi::AmanziGeometry::Point& p, double t) const
  {
    // u = a x^2 + y^2
    double x = p[0];
    double y = p[1];
    return 2 * pow(x, 2) + pow(y, 2);
  }
  double exact1(const Amanzi::AmanziGeometry::Point& p, double t) const
  {
    // v = sin^2( pi/4 * (x+y) )
    double x = p[0];
    double y = p[1];
    return pow(sin(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)), 2);
  }

  Amanzi::AmanziGeometry::Point gradient_exact0(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::AmanziGeometry::Point g(2);
    g[0] = 4 * p[0];
    g[1] = 2 * p[1];
    return g;
  }

  Amanzi::AmanziGeometry::Point gradient_exact1(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    Amanzi::AmanziGeometry::Point g(2);
    double x = p[0];
    double y = p[1];
    g[0] = (1.0L / 2.0L) * M_PI * sin(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)) *
           cos(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y));
    g[1] = (1.0L / 2.0L) * M_PI * sin(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)) *
           cos(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y));
    return g;
  }

  double source_exact0(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    double y = p[1];
    return -4 * x *
             (4 * x * pow(sin(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)), 2) +
              (1.0L / 2.0L) * M_PI * (2 * pow(x, 2) + pow(y, 2)) *
                sin(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)) *
                cos(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y))) *
             exp((2 * pow(x, 2) + pow(y, 2)) *
                 pow(sin(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)), 2)) -
           2 * y *
             (2 * y * pow(sin(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)), 2) +
              (1.0L / 2.0L) * M_PI * (2 * pow(x, 2) + pow(y, 2)) *
                sin(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)) *
                cos(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y))) *
             exp((2 * pow(x, 2) + pow(y, 2)) *
                 pow(sin(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)), 2)) -
           6 * exp((2 * pow(x, 2) + pow(y, 2)) *
                   pow(sin(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)), 2));
  }

  double source_exact1(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    double x = p[0];
    double y = p[1];
    return -1.0L / 2.0L * M_PI *
             (2 * x + (1.0L / 4.0L) * M_PI * sin(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)) *
                        cos(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y))) *
             exp(pow(x, 2) + (1.0L / 2.0L) * pow(y, 2) +
                 (1.0L / 2.0L) * pow(sin(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)), 2)) *
             sin(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)) *
             cos(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)) -
           1.0L / 2.0L * M_PI *
             (y + (1.0L / 4.0L) * M_PI * sin(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)) *
                    cos(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y))) *
             exp(pow(x, 2) + (1.0L / 2.0L) * pow(y, 2) +
                 (1.0L / 2.0L) * pow(sin(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)), 2)) *
             sin(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)) *
             cos(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)) +
           (1.0L / 4.0L) * pow(M_PI, 2) *
             exp(pow(x, 2) + (1.0L / 2.0L) * pow(y, 2) +
                 (1.0L / 2.0L) * pow(sin(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)), 2)) *
             pow(sin(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)), 2) -
           1.0L / 4.0L * pow(M_PI, 2) *
             exp(pow(x, 2) + (1.0L / 2.0L) * pow(y, 2) +
                 (1.0L / 2.0L) * pow(sin(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)), 2)) *
             pow(cos(M_PI * ((1.0L / 4.0L) * x + (1.0L / 4.0L) * y)), 2);
  }
};

#endif
