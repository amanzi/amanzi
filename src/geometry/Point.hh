/*
This is the geometry component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Rao Garimella
         Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef   __POINT_HXX
#define   __POINT_HXX

#include <iostream>
#include <vector>
#include <cmath>


namespace Amanzi {
namespace AmanziGeometry {

class Point {
 public:
  Point() {
    d = 0;
    xyz = NULL;
  }
  Point(const Point& p) {
    d = 0;
    xyz = NULL;
    init(p.dim());
    for (int i = 0; i < d; i++) xyz[i] = p.xyz[i];
  }
  explicit Point(const int N) {
    d = 0;
    xyz = NULL;
    init(N);
  }
  Point(const double& x, const double& y) {
    d = 0;
    xyz = NULL;
    init(2);
    xyz[0] = x;
    xyz[1] = y;
  }
  Point(const double& x, const double& y, const double& z) {
    d = 0;
    xyz = NULL;
    init(3);
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
  }
  ~Point() { if (xyz != NULL) delete [] xyz; }

  // main members
  void init(const int N) {
    if (!xyz) {
      xyz = new double[N];
    } else if (d < N) {  // If it was reinitialized to be a higher dim point
      delete [] xyz;
      xyz = new double[N];
    }
    d = N;
    for (int i = 0; i < d; i++) xyz[i] = 0.0;
  }

  void set(const double& val) {
    if (this == NULL) throw std::exception();
    for (int i = 0; i < d; i++) xyz[i] = val;
  }
  void set(const Point& p) {
    if (d != p.d) throw std::exception();
    for (int i = 0; i < d; i++) xyz[i] = p[i];
  }
  void set(const double *val) {
    if (val == NULL) throw std::exception();
    for (int i = 0; i < d; i++) xyz[i] = val[i];
  }
  void set(const double& x, const double& y) {
    if (d != 2) throw std::exception();
    xyz[0] = x;
    xyz[1] = y;
  }
  void set(const double& x, const double& y, const double& z) {
    if (d != 3) throw std::exception();
    xyz[0] = x;
    xyz[1] = y;
    xyz[2] = z;
  }

  int is_valid() { return (d == 2 || d == 3) ? 1 : 0; }

  // access members
  double& operator[] (const int i) { return xyz[i]; }
  const double& operator[] (const int i) const { return xyz[i]; }

  double x() const { return xyz[0]; }
  double y() const { return xyz[1]; }
  double z() const { return (this->d == 3) ? xyz[2] : 0.0; }

  int dim() const { return d; }

  // operators
  Point& operator=(const Point& p) {
    init(p.dim());
    for (int i = 0; i < d; i++) xyz[i] = p[i];
    return *this;
  }

  Point& operator+=(const Point& p) {
    for (int i = 0; i < d; i++) xyz[i] += p[i];
    return *this;
  }
  Point& operator-=(const Point& p) {
    for (int i = 0; i < d; i++) xyz[i] -= p[i];
    return *this;
  }
  Point& operator*=(const double& c) {
    for (int i = 0; i < d; i++) xyz[i] *= c;
    return *this;
  }
  Point& operator/=(const double& c) {
    for (int i = 0; i < d; i++) xyz[i] /= c;
    return *this;
  }

  friend Point operator*(const double& r, const Point& p) {
    return (p.d == 2) ? Point(r*p[0], r*p[1]) : Point(r*p[0], r*p[1], r*p[2]);
  }
  friend Point operator*(const Point& p, const double& r) { return r*p; }
  friend double operator*(const Point& p, const Point& q) {
    double s = 0.0; 
    for (int i = 0; i < p.d; i++ ) s += p[i]*q[i];
    return s;
  }

  friend Point operator/(const Point& p, const double& r) { return p * (1.0/r); }

  friend Point operator+(const Point& p, const Point& q) {
    return (p.d == 2) ? Point(p[0]+q[0], p[1]+q[1]) : Point(p[0]+q[0], p[1]+q[1], p[2]+q[2]);
  }
  friend Point operator-(const Point& p, const Point& q) {
    return (p.d == 2) ? Point(p[0]-q[0], p[1]-q[1]) : Point(p[0]-q[0], p[1]-q[1], p[2]-q[2]);
  }
  friend Point operator-(const Point& p) {
    return (p.d == 2) ? Point(-p[0], -p[1]) : Point(-p[0], -p[1], -p[2]);
  }

  friend Point operator^(const Point& p, const Point& q) {
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

  /* miscellaneous */
  inline friend double L22(const Point& p) { return p*p; }
  inline friend double norm(const Point& p) { return sqrt(p*p); }

  friend std::ostream& operator<<(std::ostream& os, const Point& p) {
    os << p.x() << " " << p.y();
    if (p.d == 3) os << " " << p.z();
    return os;
  }

 private:
  int d;
  double *xyz;
};  // end class Point

typedef std::vector<Point> Point_List;


}  // end namespace AmanziGeometry
}  // end namespace Amanzi

#endif

