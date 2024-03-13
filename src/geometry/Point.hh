/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Geometry

*/

#ifndef AMANZI_GEOMETRY_POINT_HH_
#define AMANZI_GEOMETRY_POINT_HH_

#include <iostream>
#include <vector>
#include <cmath>

#include "dbc.hh"

namespace Amanzi {
namespace AmanziGeometry {

class Point {
 public:
  Point()
  {
    d = 0;
    xyz[0] = xyz[1] = xyz[2] = 0.0;
  }
  Point(const Point& p)
  {
    d = p.d;
    std::copy(p.xyz, p.xyz + d, xyz);
  }
  explicit Point(const int N)
  {
    d = N;
    xyz[0] = xyz[1] = xyz[2] = 0.0;
  }
  Point(const double& x, const double& y)
  {
    d = 2;
    xyz[0] = x;
    xyz[1] = y;
  }
  Point(const double& x, const double& y, const double& z)
  {
    d = 3;
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
  }
  ~Point(){};

  // main members
  void set(const double& val)
  {
    AMANZI_ASSERT(d > 0);
    for (int i = 0; i < d; i++) xyz[i] = val;
  }
  void set(const Point& p)
  {
    d = p.d;
    std::copy(p.xyz, p.xyz + d, xyz);
  }
  void set(const double* val)
  {
    AMANZI_ASSERT(val);
    AMANZI_ASSERT(d > 0);
    std::copy(val, val + d, xyz);
  }
  void set(const int N, const double* val)
  {
    AMANZI_ASSERT(val);
    d = N;
    std::copy(val, val + d, xyz);
  }
  void set(const double& x, const double& y)
  {
    d = 2;
    xyz[0] = x;
    xyz[1] = y;
  }
  void set(const double& x, const double& y, const double& z)
  {
    d = 3;
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
  }

  int is_valid() { return (d == 2 || d == 3) ? 1 : 0; }

  // access members
  double& operator[](const int i) { return xyz[i]; }
  const double& operator[](const int i) const { return xyz[i]; }

  double x() const { return xyz[0]; }
  double y() const { return xyz[1]; }
  double z() const { return (d == 3) ? xyz[2] : 0.0; }

  int dim() const { return d; }

  // operators
  Point& operator=(const Point& p)
  {
    d = p.d;
    std::copy(p.xyz, p.xyz + d, xyz);
    return *this;
  }

  Point& operator+=(const Point& p)
  {
    for (int i = 0; i < d; i++) xyz[i] += p[i];
    return *this;
  }
  Point& operator-=(const Point& p)
  {
    for (int i = 0; i < d; i++) xyz[i] -= p[i];
    return *this;
  }
  Point& operator*=(const double& c)
  {
    for (int i = 0; i < d; i++) xyz[i] *= c;
    return *this;
  }
  Point& operator/=(const double& c)
  {
    for (int i = 0; i < d; i++) xyz[i] /= c;
    return *this;
  }

  friend Point operator*(const double& r, const Point& p)
  {
    return (p.d == 2) ? Point(r * p[0], r * p[1]) : Point(r * p[0], r * p[1], r * p[2]);
  }
  friend Point operator*(const Point& p, const double& r) { return r * p; }
  friend double operator*(const Point& p, const Point& q)
  {
    double s = 0.0;
    for (int i = 0; i < p.d; i++) s += p[i] * q[i];
    return s;
  }

  friend Point operator/(const Point& p, const double& r) { return p * (1.0 / r); }

  friend Point operator+(const Point& p, const Point& q)
  {
    return (p.d == 2) ? Point(p[0] + q[0], p[1] + q[1]) :
                        Point(p[0] + q[0], p[1] + q[1], p[2] + q[2]);
  }
  friend Point operator-(const Point& p, const Point& q)
  {
    return (p.d == 2) ? Point(p[0] - q[0], p[1] - q[1]) :
                        Point(p[0] - q[0], p[1] - q[1], p[2] - q[2]);
  }
  friend Point operator-(const Point& p)
  {
    return (p.d == 2) ? Point(-p[0], -p[1]) : Point(-p[0], -p[1], -p[2]);
  }

  friend Point operator^(const Point& p, const Point& q)
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
inline double
L22(const Point& p)
{
  return p * p;
}
inline double
norm(const Point& p)
{
  return sqrt(p * p);
}

inline bool
operator==(const Point& p, const Point& q)
{
  if (p.dim() != q.dim()) return false;
  for (int i = 0; i < p.dim(); ++i)
    if (p[i] != q[i]) return false;
  return true;
}

inline bool
operator!=(const Point& p, const Point& q)
{
  return !(p == q);
}

} // namespace AmanziGeometry
} // namespace Amanzi

#endif
