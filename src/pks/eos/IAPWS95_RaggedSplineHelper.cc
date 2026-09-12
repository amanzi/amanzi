/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <algorithm>
#include <iostream>
#include <limits>

#include "dbc.hh"

#include "IAPWS95_RaggedSplineHelper.hh"

namespace Amanzi {
namespace AmanziEOS {

/* ******************************************************************
* Constructs knot vector for a cubic B-spline from the coordinate lines.
* For cubic splines, the degree is p=3. A clamped knot vector repeats 
* the first and last knot p+1=4 times. This makes the spline interpolate 
* the endpoint behavior in the usual open-knot-vector sense and ensures 
* that the basis spans the full interval.
****************************************************************** */
std::vector<double>
IAPWS95_RaggedSplineHelper::MakeClampedCubicKnots(const std::vector<double>& lines)
{
  std::vector<double> knots;
  knots.reserve(lines.size() + 6);
  for (int k = 0; k < 4; ++k) knots.push_back(lines.front());
  for (int i = 1; i + 1 < lines.size(); ++i) knots.push_back(lines[i]);
  for (int k = 0; k < 4; ++k) knots.push_back(lines.back());
  return knots;
}


/* ******************************************************************
*
****************************************************************** */
IAPWS95_RaggedSplineHelper::BasisData
IAPWS95_RaggedSplineHelper::EvaluateCubicBasis(const std::vector<double>& U, double x) const
{
  const int p = 3;
  const int n = static_cast<int>(U.size()) - p - 2;

  x = std::clamp(x, U[p], U[n + 1]);
  int span;
  if (x >= U[n + 1]) {
    span = n;
  } else {
    int low = p, high = n + 1;
    span = (low + high) / 2;
    while (x < U[span] || x >= U[span + 1]) {
      if (x < U[span]) high = span;
      else low = span;
      span = (low + high) / 2;
    }
  }

  // Algorithm A2.3 from The NURBS Book: derivatives of nonzero basis functions.
  double ndu[4][4] = {};
  double left[4] = {}, right[4] = {};
  ndu[0][0] = 1.0;
  for (int j = 1; j <= p; ++j) {
    left[j] = x - U[span + 1 - j];
    right[j] = U[span + j] - x;
    double saved = 0.0;
    for (int r = 0; r < j; ++r) {
      ndu[j][r] = right[r + 1] + left[j - r];
      const double temp = ndu[r][j - 1] / ndu[j][r];
      ndu[r][j] = saved + right[r + 1] * temp;
      saved = left[j - r] * temp;
    }
    ndu[j][j] = saved;
  }

  double ders[3][4] = {};
  for (int j = 0; j <= p; ++j) ders[0][j] = ndu[j][p];

  double a[2][4] = {};
  for (int r = 0; r <= p; ++r) {
    int s1 = 0, s2 = 1;
    a[0][0] = 1.0;
    for (int k = 1; k <= 2; ++k) {
      double d = 0.0;
      const int rk = r - k;
      const int pk = p - k;
      if (r >= k) {
        a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
        d = a[s2][0] * ndu[rk][pk];
      }
      const int j1 = (rk >= -1) ? 1 : -rk;
      const int j2 = (r - 1 <= pk) ? k - 1 : p - r;
      for (int j = j1; j <= j2; ++j) {
        a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
        d += a[s2][j] * ndu[rk + j][pk];
      }
      if (r <= pk) {
        a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
        d += a[s2][k] * ndu[r][pk];
      }
      ders[k][r] = d;
      std::swap(s1, s2);
    }
  }

  int factor = p;
  for (int k = 1; k <= 2; ++k) {
    for (int j = 0; j <= p; ++j) ders[k][j] *= factor;
    factor *= (p - k);
  }

  BasisData out;
  for (int j = 0; j < 4; ++j) {
    out.index[j] = span - p + j;
    out.value[j] = ders[0][j];
    out.d1[j] = ders[1][j];
    out.d2[j] = ders[2][j];
  }
  return out;
}


/* ******************************************************************
*
****************************************************************** */
IAPWS95_RaggedSplineHelper::BasisData
IAPWS95_RaggedSplineHelper::EvaluateCachedCubicBasis(const std::vector<double>& U,
                                                     const std::vector<CubicSpanData>& cache,
                                                     const SpanLookup& lookup,
                                                     double x) const
{
  constexpr int p = 3;
  const int n = static_cast<int>(U.size()) - p - 2;

  x = std::clamp(x, U[p], U[n + 1]);

  int span = FindSpanFast(U, lookup, x);
  const CubicSpanData& S = cache[span];

  const double xi = (x - S.x0) * S.inv_h;
  double inv_h2 = S.inv_h * S.inv_h;

  BasisData out;

  for (int r = 0; r < 4; ++r) {
    const double a0 = S.poly[r][0];
    const double a1 = S.poly[r][1];
    const double a2 = S.poly[r][2];
    const double a3 = S.poly[r][3];

    out.index[r] = S.first_basis + r;

    // Horner evaluation
    out.value[r] = ((a3 * xi + a2) * xi + a1) * xi + a0;
    out.d1[r] = ((3 * a3 * xi + 2 * a2) * xi + a1) * S.inv_h;
    out.d2[r] = (6 * a3 * xi + 2 * a2) * inv_h2;
  }

  return out;
}


/* ******************************************************************
* Evaluate residual part using dimensionless input
****************************************************************** */
std::array<double, 6>
IAPWS95_RaggedSplineHelper::Evaluate(const Mesh& mesh,
                                     const std::vector<double>& coefficients,
                                     double delta, double tau) const
{
  const BasisData& Bd = EvaluateCachedCubicBasis(mesh.x_knots,
                                                 mesh.x_span_data,
                                                 mesh.x_lookup,
                                                 delta);
  const BasisData& Bt = EvaluateCachedCubicBasis(mesh.y_knots,
                                                 mesh.y_span_data,
                                                 mesh.y_lookup,
                                                 tau);
  std::array<double, 6> out{};

  for (int a = 0; a < 4; ++a) {
    int i = Bd.index[a];
    for (int c = 0; c < 4; ++c) {
      int j = Bt.index[c];
      double z = coefficients[mesh.CoefficientIndex(i, j)];
      out[0] += z * Bd.value[a] * Bt.value[c];
      out[1] += z * Bd.d1[a] * Bt.value[c];
      out[2] += z * Bd.value[a] * Bt.d1[c];
      out[3] += z * Bd.d2[a] * Bt.value[c];
      out[4] += z * Bd.d1[a] * Bt.d1[c];
      out[5] += z * Bd.value[a] * Bt.d2[c];
    }
  }
  return out;
}


std::array<double, 3>
IAPWS95_RaggedSplineHelper::EvaluateFirst(const Mesh& mesh,
                                          const std::vector<double>& coefficients,
                                          double delta, double tau) const
{
  const BasisData& Bd = EvaluateCachedCubicBasis(mesh.x_knots,
                                                 mesh.x_span_data,
                                                 mesh.x_lookup,
                                                 delta);
  const BasisData& Bt = EvaluateCachedCubicBasis(mesh.y_knots,
                                                 mesh.y_span_data,
                                                 mesh.y_lookup,
                                                 tau);
  std::array<double, 3> out{};

  for (int a = 0; a < 4; ++a) {
    int i = Bd.index[a];
    for (int c = 0; c < 4; ++c) {
      int j = Bt.index[c];
      double z = coefficients[mesh.CoefficientIndex(i, j)];
      out[0] += z * Bd.value[a] * Bt.value[c];
      out[1] += z * Bd.d1[a] * Bt.value[c];
      out[2] += z * Bd.value[a] * Bt.d1[c];
    }
  }
  return out;
}


/* ******************************************************************
* Compute cache data for all spans
****************************************************************** */
std::vector<IAPWS95_RaggedSplineHelper::CubicSpanData>
IAPWS95_RaggedSplineHelper::BuildCubicSpanCache(const std::vector<double>& U)
{
  constexpr int p = 3;
  const int n = static_cast<int>(U.size()) - p - 2;
  AMANZI_ASSERT(n >= p);  // Invalid cubic knot vector.

  std::vector<CubicSpanData> cache(n + 1);

  constexpr double xa = 0.25;
  constexpr double xb = 0.75;

  for (int span = p; span <= n; ++span) {
    double x0 = U[span];
    double x1 = U[span + 1];
    // Repeated knots have zero-width spans.
    if (!(x1 > x0)) continue;

    double h = x1 - x0;
    CubicSpanData& S = cache[span];

    S.first_basis = span - p;
    S.x0 = x0;
    S.inv_h = 1.0 / h;

    // Evaluate the existing reference implementation at xi = 1/4 and xi = 3/4.
    const BasisData Ba = EvaluateCubicBasis(U, x0 + xa * h);
    const BasisData Bb = EvaluateCubicBasis(U, x0 + xb * h);

    for (int r = 0; r < 4; ++r) {
      double fa = Ba.value[r];
      double fb = Bb.value[r];

      // convert derivatives wrt to x into derivatives wrt to xi.
      double da = h * Ba.d1[r];
      double db = h * Bb.d1[r];

      S.poly[r][0] = fb - (9.0 / 16.0) * da - (3.0 / 16.0) * db;
      S.poly[r][1] = 9.0 * (fa - fb) + (15.0 / 4.0) * da + (7.0 / 4.0) * db;
      S.poly[r][2] = 24.0 * (fb - fa) - 7.0 * da - 5.0 * db;
      S.poly[r][3] = 16.0 * (fa - fb) + 4.0 * (da + db);
    }
  }

  return cache;
}


/* ******************************************************************
* Minimum dT
****************************************************************** */
double
IAPWS95_RaggedSplineHelper::MinSpacing(const std::vector<double>& x) const
{
  double h = std::numeric_limits<double>::max();
  for (int i = 0; i + 1 < x.size(); ++i) {
    h = std::min(h, x[i + 1] - x[i]);
  }
  return h;
}


/* ******************************************************************
*
****************************************************************** */
int
IAPWS95_RaggedSplineHelper::LowerCell(const std::vector<double>& x, double value)
{
  AMANZI_ASSERT(x.size() >= 2); // Coordinate line has fewer than two points.
  if (value <= x.front()) return 0;
  if (value >= x.back()) return x.size() - 2;
  return std::distance(x.begin(), std::upper_bound(x.begin(), x.end(), value)) - 1;
}


/* ******************************************************************
* Build a uniform auxiliary lookup table for a non-uniform coordinate
* mesh x[0] < x[1] < ... < x[n].
*
* The entry lookup.cell[b] contains the mesh cell containing the
* LEFT EDGE of bin b.
*
* Therefore, for any point inside bin b, the stored cell can only be
* equal to or to the left of the actual cell.
****************************************************************** */
IAPWS95_RaggedSplineHelper::SpanLookup
IAPWS95_RaggedSplineHelper::MakeSpanLookupLeft(const std::vector<double>& U, int nbins)
{
  constexpr int p = 3;
  const int n = static_cast<int>(U.size()) - p - 2;

  SpanLookup lookup;

  lookup.xmin = U[p];
  lookup.xmax = U[n + 1];
  lookup.nbins = nbins;
  lookup.scale = static_cast<double>(nbins) / (lookup.xmax - lookup.xmin);

  lookup.span.resize(nbins);

  int span = p;
  for (int b = 0; b < nbins; ++b) {
    // Left edge of auxiliary bin.
    double xb = lookup.xmin + static_cast<double>(b) / lookup.scale;

    // Find spline span satisfying  U[span] <= xb < U[span+1].
    // Span starts at p, so the repeated knots U[0],...,U[p] do not enter the search.
    while (span < n && xb >= U[span + 1]) span++;
    lookup.span[b] = span;
  }

  return lookup;
}


/* ******************************************************************
* Fast lookup of mesh cell containing value: x[i] <= value < x[i+1].
* At the right endpoint return the last cell.
* Since lookup.cell[b] corresponds to the LEFT edge of the bin,
* correction is required only in the increasing direction.
****************************************************************** */
int
IAPWS95_RaggedSplineHelper::FindSpanFast(const std::vector<double>& U,
                                         const SpanLookup& lookup,
                                         double x) const
{
  constexpr int p = 3;
  const int n = static_cast<int>(U.size()) - p - 2;

  if (x <= lookup.xmin) return p;
  if (x >= lookup.xmax) return n;

  int b = static_cast<int>((x - lookup.xmin) * lookup.scale);
  int span = lookup.span[b];

  // Because the table is based on the LEFT bin edge,
  // its span cannot be to the right of the true span.
  while (span < n && x >= U[span + 1]) span++;
  return span;
}

}  // namespace AmanziEOS
}  // namespace Amanzi
