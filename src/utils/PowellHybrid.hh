/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

#ifndef UTILS_POWELL_HYBRID_HH_
#define UTILS_POWELL_HYBRID_HH_

#include <vector>
#include <cmath>

#include "dbc.hh"


/* ******************************************************************
* Powell hybrid method for solving a system of equtions.
*
* Two temination criteria are used, tolerance and the maximum number
* of iterations.
*
* Upon return:
*   *itr == -1 indicates that method did not converge
*   otherwise, *itr is the number of performed iterations.
****************************************************************** */
template <class Function, class Vector>
Vector PowellHybrid(Vector& x, const Function& f, int* itrs, double tol = 1e-10);

template<class Matrix, class Vector>
Vector DoglegStep(const Vector& b, const Matrix& J, double delta);

template<class Function, class Matrix, class Vector>
void FiniteDifferenceJacobian(const Vector& x, const Vector& Fx, Function& Ffun, Matrix& J);

template<class Matrix, class Vector>
Vector MatrixVector(const Matrix& A, const Vector& x);

template<class Matrix, class Vector>
Vector GaussSolve(Matrix A, Vector b);

template<class Matrix>
Matrix Transpose(const Matrix& A);


/* ******************************************************************
* Powell hybrid method.
****************************************************************** */
template <class Function, class Vector>
Vector PowellHybrid(Vector& x, const Function& fun, int* itrs, double tol)
{
  using Matrix = std::vector<Vector>;

  int n = x.size();
  int maxit = *itrs;
  double delta = 1.0;

  Vector r(x), rnew(x), rpred(x);
  Matrix J(n, r);

  for (int k = 0; k < maxit; ++k) {
    r = fun(x);
    double rnorm = norm(r);

    if (rnorm < tol) {
      *itrs = k;
      return x;
    }

    FiniteDifferenceJacobian(x, r, fun, J);
    Vector s = DoglegStep(r, J, delta);

    Vector x_new = x + s;
    rnew = fun(x_new);
    double rnorm_new = norm(rnew);
    double ared = rnorm * rnorm - rnorm_new * rnorm_new;

    rpred = r + MatrixVector(J, s);
    double rnorm_pred = norm(rpred);

    double pred = rnorm * rnorm - rnorm_pred * rnorm_pred;
    double rho = ared / pred;

    if (rho > 0.1) {
      x = x_new;
      if (rho > 0.75) delta = std::max(delta, 2.0 * norm(s));
    } else {
      delta *= 0.5;
    }
  }

  *itrs = -1;
  return x;
}


/* ******************************************************************
* Dogleg strategy to chose a step inside a trust region.
****************************************************************** */
template<class Matrix, class Vector>
Vector DoglegStep(const Vector& f, const Matrix& J, double delta)
{
  int n = f.size();

  // Newton (sN) step
  Vector sN = GaussSolve(J, (-1.0) * f);
  if (norm(sN) <= delta) return sN;

  // steepest descent (-g) step
  Matrix JT = Transpose(J);
  Vector g = MatrixVector(JT, f);

  double gnorm = norm(g);
  if (gnorm < 1e-14) return (-delta / gnorm) * g;

  double tmp = norm(MatrixVector(J, g));
  double alpha = (gnorm * gnorm) / (tmp * tmp);

  Vector sU = (-alpha) * g;

  if (norm(sU) >= delta) return (delta / norm(sU)) * sU;

  // interpolate between sU and sN
  Vector d = sN - sU;

  double a(0.0), b(0.0), c(0.0);
  for (int i = 0; i < n; ++i) {
    a += d[i] * d[i];
    b += 2.0 * sU[i] * d[i];
    c += sU[i] * sU[i];
  }
  c -= delta * delta;

  double tau = (-b + std::sqrt(b * b - 4 * a * c)) / (2 * a);
  return sU + tau * d;
}


/* ******************************************************************
* Finite differece approximation of Jacobian
****************************************************************** */
template<class Function, class Matrix, class Vector>
void FiniteDifferenceJacobian( const Vector& x, const Vector& fx, Function& f, Matrix& J)
{
  const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
  int n = x.size();

  Vector xh(x), fxh(x);

  for (int j = 0; j < n; ++j) {
    double h = eps * std::max(std::abs(x[j]), 1.0);
    xh[j] += h;

    Vector fxh = f(xh);

    for (size_t i = 0; i < n; ++i)
      J[i][j] = (fxh[i] - fx[i]) / h;

    xh[j] = x[j];
  }
}


/* ******************************************************************
* Supporting function: matrix-vector product
****************************************************************** */
template<class Matrix, class Vector>
Vector MatrixVector(const Matrix& A, const Vector& x)
{
  Vector r(x.size());
  for (size_t i = 0; i < A.size(); ++i) {
    r[i] = 0.0;
    for (size_t j = 0; j < x.size(); ++j)
      r[i] += A[i][j] * x[j];
  }
  return r;
}


/* ******************************************************************
* Supporting function: Gauss ellimination
****************************************************************** */
template<class Matrix, class Vector>
Vector GaussSolve(Matrix A, Vector b)
{
  const int n = b.size();

  for (size_t k = 0; k < n; ++k) {
    double pivot = A[k][k];
    AMANZI_ASSERT(std::fabs(pivot) > 1e-14);

    for (int j = k; j < n; ++j) A[k][j] /= pivot;
    b[k] /= pivot;

    for (int i = k + 1; i < n; ++i) {
      double f = A[i][k];
      for (size_t j = k; j < n; ++j)
        A[i][j] -= f * A[k][j];
      b[i] -= f * b[k];
    }
  }

  Vector x(n);
  for (int i = n - 1; i >= 0; --i) {
    x[i] = b[i];
    for (size_t j = i + 1; j < n; ++j)
      x[i] -= A[i][j] * x[j];
  }
  return x;
}


/* ******************************************************************
* Supporting function: matrix transpose
****************************************************************** */
template<class Matrix>
Matrix Transpose(const Matrix& A)
{
  Matrix AT(A);
  for (size_t i = 0; i < A.size(); ++i)
    for (size_t j = 0; j < A[0].size(); ++j)
      AT[j][i] = A[i][j];
  return AT;
}

#endif
