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
****************************************************************** */
template <class F>
double
computeRootBrent(const F& f, double a, double b, double tol, int* itr)
{
  int itr_max(*itr);
  double c, d, s, fa, fb, fc, fs, ftol;
  bool flag(true);

  fa = f(a);
  fb = f(b);
  if (fa * fb >= 0.0) {
    *itr = -1;
    return 0.0;
  }
  ftol = (1.0 + std::fabs(fa) + std::fabs(fb)) * tol;

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

  return 0.0; // default value
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

} // namespace Utils
} // namespace Amanzi

#endif
