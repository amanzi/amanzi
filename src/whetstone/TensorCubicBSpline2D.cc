/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  A faster implementation of not-a-knot 2D tensor-product cubic spline.
*/

#include "TensorCubicBSpline2D.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* values[j][i] = f(x[i], y[j])
****************************************************************** */
int
TensorCubicBSpline2D::build(const std::vector<std::vector<double>>& values)
{
  nx_ = x_.size();
  ny_ = y_.size();

  if (nx_ < 4 || ny_ < 4) {
    return -1; // Cubic not-a-knot spline requires nx, ny >= 4.
  }

  if (checkIncreasing(x_) != 0) return -2; // grid must be strictly increasing.
  if (checkIncreasing(y_) != 0) return -2; 

  if (values.size() != ny_) return -3;  // values must have ny rows
  for (int j = 0; j < ny_; ++j) {
    if (values[j].size() != nx_) return -4;  // Each row of values must have nx entries
  }

  tx_ = makeNotAKnotKnots_(x_);
  ty_ = makeNotAKnotKnots_(y_);

  // F is stored column-major as nx by ny:
  // F(i,j) = values[j][i]
  std::vector<double> F(nx_ * ny_);
  for (int j = 0; j < ny_; ++j) {
    for (int i = 0; i < nx_; ++i) {
      F[i + nx_ * j] = values[j][i];
    }
  }

  // Solve A_x * Y = F, Y is nx by ny
  int info = 0;
  char trans = 'N';
  std::vector<int> ipiv(std::max(nx_, ny_));
  auto Ax = buildCollocationMatrix_(x_, tx_);

  DGETRF_F77(&nx_, &nx_, Ax.data(), &nx_, ipiv.data(), &info);
  if (info < 0) return -5;  // dgetrf: illegal argument
  if (info > 0) return -6;  // dgetrf: singular matrix

  std::vector<double> Y = F;
  DGETRS_F77(&trans, &nx_, &ny_, Ax.data(), &nx_, ipiv.data(), Y.data(), &nx_, &info);
  if (info != 0) return -7;  // dgetrs failed

  // Solve A_y * C^T = Y^T, YT is ny by nx
  auto Ay = buildCollocationMatrix_(y_, ty_);

  std::vector<double> YT(ny_ * nx_);
  for (int j = 0; j < ny_; ++j) {
    for (int i = 0; i < nx_; ++i) {
      YT[j + ny_ * i] = Y[i + nx_ * j];
    }
  }

  DGETRF_F77(&ny_, &ny_, Ay.data(), &ny_, ipiv.data(), &info);
  if (info < 0) return -5;  // dgetrf: illegal argument
  if (info > 0) return -6;  // dgetrf: singular matrix

  DGETRS_F77(&trans, &ny_, &nx_, Ay.data(), &ny_, ipiv.data(), YT.data(), &ny_, &info);
  if (info != 0) return -7;  // dgetrs failed

  // Store coefficients C(j,i) in row-major coeff_[j * nx_ + i].
  coeff_.assign(nx_ * ny_, 0.0);
  for (int i = 0; i < nx_; ++i) {
    for (int j = 0; j < ny_; ++j) {
      coeff_[j * nx_ + i] = YT[j + ny_ * i];
    }
  }

  return 0;
}


/* ******************************************************************
* values[j][i] = f(x[i], y[j])
****************************************************************** */
std::array<double, 6>
TensorCubicBSpline2D::evaluate(double x, double y) const
{
  constexpr int p = 3;

  // Clamp evaluation to interpolation domain.
  x = std::min(std::max(x, x_.front()), x_.back());
  y = std::min(std::max(y, y_.front()), y_.back());

  int span_x = findSpan_(nx_, p, x, tx_);
  int span_y = findSpan_(ny_, p, y, ty_);

  auto bx = basisDers_(span_x, x, p, 2, tx_);
  auto by = basisDers_(span_y, y, p, 2, ty_);

  double S(0.0), Sx(0.0), Sy(0.0), Sxx(0.0), Sxy(0.0), Syy(0.0);

  for (int i = 0; i <= p; ++i) {
    int ix = span_x - p + i;

    double Bx = bx[0][i];
    double dBx = bx[1][i];
    double ddBx = bx[2][i];

    for (int j = 0; j <= p; ++j) {
      int iy = span_y - p + j;
      double c = coeff_[iy * nx_ + ix];

      double By = by[0][j];
      double dBy = by[1][j];
      double ddBy = by[2][j];

      S += c * Bx * By;
      Sx += c * dBx * By;
      Sy += c * Bx * dBy;
      Sxx += c * ddBx * By;
      Sxy += c * dBx * dBy;
      Syy += c * Bx * ddBy;
    }
  }

  return {S, Sx, Sy, Sxx, Sxy, Syy};
}


/* ******************************************************************
* Cubic not-a-knot knot vector for interpolation points x.
****************************************************************** */
std::vector<double>
TensorCubicBSpline2D::makeNotAKnotKnots_(const std::vector<double>& x)
{
  int n = x.size();
  constexpr int p = 3;

  std::vector<double> t(n + p + 1);

  t[0] = x.front();
  t[1] = x.front();
  t[2] = x.front();
  t[3] = x.front();

  int k = 4;
  for (int i = 2; i <= n - 3; ++i) t[k++] = x[i];

  t[k++] = x.back();
  t[k++] = x.back();
  t[k++] = x.back();
  t[k++] = x.back();

  return t;
}


/* ******************************************************************
* Build dense collocation matrix A(i,j) = B_j(x_i).
****************************************************************** */
std::vector<double>
TensorCubicBSpline2D::buildCollocationMatrix_(const std::vector<double>& x,
                                              const std::vector<double>& t)
{
  int n = x.size();
  constexpr int p = 3;

  std::vector<double> A(n * n, 0.0);

  for (int i = 0; i < n; ++i) {
    int span = findSpan_(n, p, x[i], t);
    auto ders = basisDers_(span, x[i], p, 0, t);

    for (int a = 0; a <= p; ++a) {
      int col = span - p + a;
      if (col >= 0 && col < n) {
        A[i + n * col] = ders[0][a]; // column-major
      }
    }
  }

  return A;
}


/* ******************************************************************
* Find span for n basis functions of degree p.
****************************************************************** */
int
TensorCubicBSpline2D::findSpan_(int nBasis, int p, double u, const std::vector<double>& U) const
{
  const int n = nBasis - 1;

  if (u >= U[n + 1]) return n;
  if (u <= U[p]) return p;

  int low = p;
  int high = n + 1;
  int mid = (low + high) / 2;

  while (u < U[mid] || u >= U[mid + 1]) {
    if (u < U[mid])  high = mid;
    else low = mid;

    mid = (low + high) / 2;
  }

  return mid;
}


/* ******************************************************************
* DersBasisFuns from The NURBS Book.
****************************************************************** */
std::array<std::array<double, 4>, 3>
TensorCubicBSpline2D::basisDers_(int span, double u, int p, int nDeriv,
                                 const std::vector<double>& U) const
{
  std::array<std::array<double, 4>, 3> ders{};
  for (auto& row : ders) row.fill(0.0);

  double ndu[4][4] = {};
  double left[4] = {};
  double right[4] = {};

  ndu[0][0] = 1.0;

  for (int j = 1; j <= p; ++j) {
    left[j] = u - U[span + 1 - j];
    right[j] = U[span + j] - u;

    double saved = 0.0;

    for (int r = 0; r < j; ++r) {
      ndu[j][r] = right[r + 1] + left[j - r];

      double temp = 0.0;
      if (ndu[j][r] != 0.0) {
        temp = ndu[r][j - 1] / ndu[j][r];
      }

      ndu[r][j] = saved + right[r + 1] * temp;
      saved = left[j - r] * temp;
    }

    ndu[j][j] = saved;
  }

  for (int j = 0; j <= p; ++j) {
    ders[0][j] = ndu[j][p];
  }

  const int du = std::min(nDeriv, p);
  double a[2][4] = {};

  for (int r = 0; r <= p; ++r) {
    int s1 = 0;
    int s2 = 1;

    a[0][0] = 1.0;

    for (int k = 1; k <= du; ++k) {
      double d = 0.0;
      int rk = r - k;
      int pk = p - k;

      int j1 = 0;
      int j2 = 0;

      if (r >= k) {
        if (ndu[pk + 1][rk] != 0.0) {
          a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
        } else {
          a[s2][0] = 0.0;
        }
        d = a[s2][0] * ndu[rk][pk];
      }

      if (rk >= -1) j1 = 1;
      else j1 = -rk;

      if (r - 1 <= pk) j2 = k - 1;
      else j2 = p - r;

      for (int j = j1; j <= j2; ++j) {
        int idx = rk + j;
        if (ndu[pk + 1][idx] != 0.0) {
          a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][idx];
        } else {
          a[s2][j] = 0.0;
        }
        d += a[s2][j] * ndu[idx][pk];
      }

      if (r <= pk) {
        if (ndu[pk + 1][r] != 0.0) {
          a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
        } else {
          a[s2][k] = 0.0;
        }
        d += a[s2][k] * ndu[r][pk];
      }

      ders[k][r] = d;
      std::swap(s1, s2);
    }
  }

  // Multiply by the correct factors:
  // first derivative: p
  // second derivative: p * (p - 1)
  int factor = p;
  for (int k = 1; k <= du; ++k) {
    for (int j = 0; j <= p; ++j) {
       ders[k][j] *= factor;
    }
    factor *= (p - k);
  }

  return ders;
}

}  // namespace WhetStone
}  // namespace Amanzi
