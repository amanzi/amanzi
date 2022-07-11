/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Rao Garimella
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_GEOMETRY_POINT_HH_
#define AMANZI_GEOMETRY_POINT_HH_

#include <iostream>
#include <vector>
#include <cmath>

#include "Kokkos_Core.hpp"
#include <Kokkos_ArithTraits.hpp>
#include "dbc.hh"

#include <cassert>

namespace Amanzi {
namespace AmanziGeometry {

class Point {
 public:
  KOKKOS_INLINE_FUNCTION Point()
  {
    d = 0;
    xyz[0] = xyz[1] = xyz[2] = 0.0;
  }
  KOKKOS_INLINE_FUNCTION Point(const Point& p)
  {
    d = p.d;
    xyz[0] = p.xyz[0];
    xyz[1] = p.xyz[1];
    xyz[2] = p.xyz[2];
  }
  KOKKOS_INLINE_FUNCTION explicit Point(const int N)
  {
    d = N;
    xyz[0] = xyz[1] = xyz[2] = 0.0;
  }
  KOKKOS_INLINE_FUNCTION Point(const double& x, const double& y)
  {
    d = 2;
    xyz[0] = x;
    xyz[1] = y;
  }
  KOKKOS_INLINE_FUNCTION
  Point(const double& x, const double& y, const double& z)
  {
    d = 3;
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
  }
  KOKKOS_INLINE_FUNCTION ~Point(){};

  // main members
  KOKKOS_INLINE_FUNCTION void set(const double& val)
  {
    assert(d > 0);
    for (int i = 0; i < d; i++) xyz[i] = val;
  }
  KOKKOS_INLINE_FUNCTION void set(const Point& p)
  {
    d = p.d;
    xyz[0] = p.xyz[0];
    xyz[1] = p.xyz[1];
    xyz[2] = p.xyz[2];
  }
  KOKKOS_INLINE_FUNCTION void set(const double* val)
  {
    assert(val);
    assert(d > 0);
    xyz[0] = val[0];
    if (d > 1) xyz[1] = val[1];
    if (d > 2) xyz[2] = val[2];
  }
  KOKKOS_INLINE_FUNCTION void set(const int N, const double* val)
  {
    assert(val);
    d = N;
    set(val);
  }
  KOKKOS_INLINE_FUNCTION void set(const double& x, const double& y)
  {
    d = 2;
    xyz[0] = x;
    xyz[1] = y;
  }
  KOKKOS_INLINE_FUNCTION void
  set(const double& x, const double& y, const double& z)
  {
    d = 3;
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
  }

  KOKKOS_INLINE_FUNCTION int is_valid() { return (d == 2 || d == 3) ? 1 : 0; }

  // access members
  KOKKOS_INLINE_FUNCTION double& operator[](const int i) { return xyz[i]; }
  KOKKOS_INLINE_FUNCTION const double& operator[](const int i) const
  {
    return xyz[i];
  }

  KOKKOS_INLINE_FUNCTION double x() const { return xyz[0]; }
  KOKKOS_INLINE_FUNCTION double y() const { return xyz[1]; }
  KOKKOS_INLINE_FUNCTION double z() const { return (d == 3) ? xyz[2] : 0.0; }

  KOKKOS_INLINE_FUNCTION int dim() const { return d; }

  KOKKOS_FORCEINLINE_FUNCTION Point nan()
  {
    return Point(Kokkos::ArithTraits<double>::nan(),
                 Kokkos::ArithTraits<double>::nan(),
                 Kokkos::ArithTraits<double>::nan());
  }


  // operators
  KOKKOS_INLINE_FUNCTION Point& operator=(const Point& p)
  {
    d = p.d;
    this->xyz[0] = p.xyz[0];
    this->xyz[1] = p.xyz[1];
    this->xyz[2] = p.xyz[2];
    return *this;
  }

  KOKKOS_INLINE_FUNCTION Point& operator+=(const Point& p)
  {
    for (int i = 0; i < d; i++) xyz[i] += p[i];
    return *this;
  }
  KOKKOS_INLINE_FUNCTION Point& operator-=(const Point& p)
  {
    for (int i = 0; i < d; i++) xyz[i] -= p[i];
    return *this;
  }
  KOKKOS_INLINE_FUNCTION Point& operator*=(const double& c)
  {
    for (int i = 0; i < d; i++) xyz[i] *= c;
    return *this;
  }
  KOKKOS_INLINE_FUNCTION Point& operator/=(const double& c)
  {
    for (int i = 0; i < d; i++) xyz[i] /= c;
    return *this;
  }

  friend KOKKOS_INLINE_FUNCTION Point operator*(const double& r, const Point& p)
  {
    return (p.d == 2) ? Point(r * p[0], r * p[1]) :
                        Point(r * p[0], r * p[1], r * p[2]);
  }
  friend KOKKOS_INLINE_FUNCTION Point operator*(const Point& p, const double& r)
  {
    return r * p;
  }
  friend KOKKOS_INLINE_FUNCTION double operator*(const Point& p, const Point& q)
  {
    double s = 0.0;
    for (int i = 0; i < p.d; i++) s += p[i] * q[i];
    return s;
  }

  friend KOKKOS_INLINE_FUNCTION Point operator/(const Point& p, const double& r)
  {
    return p * (1.0 / r);
  }

  friend KOKKOS_INLINE_FUNCTION Point operator+(const Point& p, const Point& q)
  {
    return (p.d == 2) ? Point(p[0] + q[0], p[1] + q[1]) :
                        Point(p[0] + q[0], p[1] + q[1], p[2] + q[2]);
  }
  friend KOKKOS_INLINE_FUNCTION Point operator-(const Point& p, const Point& q)
  {
    return (p.d == 2) ? Point(p[0] - q[0], p[1] - q[1]) :
                        Point(p[0] - q[0], p[1] - q[1], p[2] - q[2]);
  }
  friend KOKKOS_INLINE_FUNCTION Point operator-(const Point& p)
  {
    return (p.d == 2) ? Point(-p[0], -p[1]) : Point(-p[0], -p[1], -p[2]);
  }

  friend KOKKOS_INLINE_FUNCTION Point operator^(const Point& p, const Point& q)
  {
    Point pq(p.d);
    if (p.d == 2) {
      pq[0] = p[0] * q[1] - q[0] * p[1];
    } else if (p.d == 3) {
      pq[0] = p[1] * q[2] - p[2] * q[1];
      pq[1] = p[2] * q[0] - p[0] * q[2];
      pq[2] = p[0] * q[1] - p[1] * q[0];
    }
    return pq;
  }

  friend std::ostream& operator<<(std::ostream& os, const Point& p)
  {
    os << p.x() << " " << p.y();
    if (p.d == 3) os << " " << p.z();
    return os;
  }

 private:
  int d;
  double xyz[3];
}; // class Point


// Miscellaneous non-member functions
KOKKOS_INLINE_FUNCTION double
L22(const Point& p)
{
  return p * p;
}

KOKKOS_INLINE_FUNCTION double
norm(const Point& p)
{
  return sqrt(p * p);
}

KOKKOS_INLINE_FUNCTION bool
operator==(const Point& p, const Point& q)
{
  if (p.dim() != q.dim()) return false;
  for (int i = 0; i < p.dim(); ++i)
    if (p[i] != q[i]) return false;
  return true;
}

KOKKOS_INLINE_FUNCTION bool
operator!=(const Point& p, const Point& q)
{
  return !(p == q);
}

} // namespace AmanziGeometry
} // namespace Amanzi

#endif
