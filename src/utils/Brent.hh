/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

#ifndef UTILS_BRENT_HH_
#define UTILS_BRENT_HH_

#include <algorithm>
#include <cmath>
#include <iostream>
#include "dbc.hh"

namespace Amanzi {
namespace Utils {

/* ******************************************************************
* Brent's method of root finding
*
* Note that a and b must bracket the root, i.e. f(a) * f(b) < 0
*
* Two temination criteria are used. First, the bracket size is less 
* than tol. Second, function values at bracket ends are ess then ftol.
* When ftol is missing or negative, it is approximated inside as 
*
*   ftol = (1 + |f(a)| + |f(b)|) * tol.
*
* Upon return:
*   *itr == -1 indicates that [a,b] did not bracket a root
*   *itr == *itr + 1 indicates that the method did not converge.
*   otherwise, *itr is the number of iterations to convergence.
****************************************************************** */
template <class F>
double
findRootBrent(const F& f, double a, double b, double tol, int* itr, double ftol = -1.0)
{
  AMANZI_ASSERT(*itr > 0);
  int itr_max(*itr);
  double c, d, s, fa, fb, fc, fs;
  bool flag(true);

  fa = f(a);
  fb = f(b);
  if (fa * fb >= 0.0) {
    *itr = -1;
    return 0.0;
  }
  if (ftol < 0.0) ftol = (1.0 + std::fabs(fa) + std::fabs(fb)) * tol;

  if (std::fabs(fa) < std::fabs(fb)) {
    std::swap(a, b);
    std::swap(fa, fb);
  }

  c = a;
  fc = fa;

  *itr = 0;
  while (*itr < itr_max) {
    (*itr)++;
    if (fa != fc && fb != fc) {
      s = a * fb * fc / ((fa - fb) * (fa - fc)) + b * fa * fc / ((fb - fa) * (fb - fc)) +
          c * fa * fb / ((fc - fa) * (fc - fb));
    } else {
      s = (fa * b - fb * a) / (fa - fb);
    }

    if ((s - (3 * a + b) / 4) * (s - b) >= 0.0 ||
        (flag && std::fabs(s - b) >= std::fabs(b - c) / 2) ||
        (!flag && std::fabs(s - b) >= std::fabs(c - d) / 2) || (flag && std::fabs(b - c) < tol) ||
        (!flag && std::fabs(c - d) < tol)) {
      s = (a + b) / 2;
      flag = true;
    } else {
      flag = false;
    }

    fs = f(s);
    d = c;
    c = b;
    fc = fb;

    if (fa * fs < 0.0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }

    if (std::fabs(fa) < std::fabs(fb)) {
      std::swap(a, b);
      std::swap(fa, fb);
    }
    if (std::fabs(fb) <= ftol) return b;
    if (std::fabs(fs) <= ftol) return s;
    if (std::fabs(b - a) < tol) return s;
  }

  (*itr)++; // indicate nonconvergence
  return s;
}


/* ******************************************************************
* Bracket a root for use in findRootBrent()
*
* Return <a,b> where f(a) * f(b) is negative, given a starting point
****************************************************************** */
template <class F>
std::pair<double, double>
bracketRoot(const F& f, double start, double delta, int* itrs)
{
  AMANZI_ASSERT(delta > 0.);
  AMANZI_ASSERT(itrs != nullptr);
  AMANZI_ASSERT(*itrs >= 0);
  double a = start;
  double b = start + delta;

  double fa = f(a);
  double fb = f(b);

  int max_itrs(*itrs);
  *itrs = 0;
  while (*itrs < max_itrs) {
    if (fa == 0.)
      return std::make_pair(a, a);
    else if (fb == 0.)
      return std::make_pair(b, b);
    else if (fa * fb < 0.)
      return std::make_pair(a, b);
    else if (std::fabs(fa) > std::fabs(fb)) {
      // root to the right of b
      std::swap(a, b);
      std::swap(fa, fb);
      b = a + delta;
      fb = f(b);
    } else {
      // root to the left of a
      std::swap(a, b);
      std::swap(fa, fb);
      a = b - delta;
      fa = f(a);
    }
    ++(*itrs);
  }
  return std::make_pair(a, b);
}


/* ******************************************************************
* A local minimum of a function f(x) in an interval [a, b].
* Follows closely to https://www.youtube.com/watch?v=BQm7uTYC0sg
****************************************************************** */
template <class F>
double
findMinimumBrent(const F& f, double a, double b, double tol, int* itr)
{
  int itr_max(*itr);
  double ratio(0.381966011250105); // golden ratio
  double x, w, v, u, e(0.0), c, xtol, xtol2, d, r, q, p;
  double fx, fw, fv, fu;

  v = w = x = b; // a + ratio * (b - a);
  fv = fw = fx = f(x);

  *itr = 0;
  while (*itr < itr_max) {
    xtol = tol * std::fabs(x) + tol;
    xtol2 = 2 * xtol;

    c = (a + b) / 2;  // called m in the video
    if (std::fabs(x - c) <= xtol2 - (b - a) / 2) return c;

    (*itr)++;

    // parabolic fit to f(x), f(v), f(w) 
    // SPI behavior is poor if any of the three conditions holds 
    //   q = 0
    //   u is ourside of [a, b]
    //   current p/q > 1/2 or previous p/q 
    p = q = r = 0.0;

    if (std::fabs(e) > xtol) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0) p = -p;
      q = std::fabs(q);
      r = e;
      e = d;
    }

    // parabolic interpolation
    if (std::fabs(p) < std::fabs(0.5 * q * r) && q * (a - x) < p && p < q * (b - x)) {
      d = p / q;
      u = x + d;

      // minimum estimate is too close to the end points
      if ((u - a) < xtol2 || (b - u) < xtol2) 
        d = (x < c) ? xtol : -xtol;

    // golden-section interpolation
    } else {
      e = (x < c) ? b - x : a - x;
      d = ratio * e;
    }

    // minumum estimate is too close to end points
    if (xtol <= std::fabs(d)) {
      u = x + d;
    } else if (d > 0.0) {
      u = x + xtol;
    } else {
      u = x - xtol;
    }

    fu = f(u);

    // update variables
    if (fu <= fx) {
      if (u < x)
        b = x;
      else
        a = x;

      v = w;
      w = x;
      x = u;
      fv = fw;
      fw = fx;
      fx = fu;
    } else {
      if (u < x)
        a = u;
      else
        b = u;

      if (fu <= fw || w == x) {
        v = w;
        w = u;
        fv = fw;
        fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u;
        fv = fu;
      }
    }
  }

  (*itr)++; // indicate nonconvergence
  return x;
}

} // namespace Utils
} // namespace Amanzi

#endif
