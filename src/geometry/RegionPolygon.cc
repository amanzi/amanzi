/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  A closed polygonal segment of a plane.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
*/

#include "dbc.hh"
#include "errors.hh"

#include "Point.hh"
#include "RegionPolygon.hh"

namespace Amanzi {
namespace AmanziGeometry {

//
// Polygon:: constructor
// -------------------------------------------------------------
RegionPolygon::RegionPolygon(const std::string& name,
                             const int id,
                             const std::vector<Point>& points,
                             const LifeCycleType lifecycle)
  : Region(name, id, true, RegionType::POLYGON, points[0].dim() - 1, points[0].dim(), lifecycle),
    normal_(points[0].dim()),
    elim_dir_(0)
{
  for (std::vector<Point>::const_iterator itr = points.begin(); itr != points.end(); ++itr) {
    points_.push_back(*itr);
  }

  Init_();
}


void
RegionPolygon::Init_()
{
  if (size() < get_manifold_dimension()) {
    Errors::Message mesg;
    mesg << "Polygons of dimension " << (int)get_manifold_dimension()
         << " need to be specified by at least " << (int)get_manifold_dimension() << " points";
    Exceptions::amanzi_throw(mesg);
  }

  int space_dimension = points_[0].dim();
  if (space_dimension == 2) {
    Point vec = points_[1] - points_[0];
    vec /= norm(vec);
    normal_.set(vec[1], -vec[0]);

    elim_dir_ = (vec[0] < vec[1]) ? 0 : 1;

  } else if (space_dimension == 3) {
    Point vec0 = points_[2] - points_[1];
    Point vec1 = points_[0] - points_[1];
    normal_ = vec0 ^ vec1;
    normal_ /= norm(normal_);

    for (int i = 3; i != points_.size(); ++i) {
      vec0 = points_[(i + 1) % size()] - points_[i];
      vec1 = points_[(i - 1 + size()) % size()] - points_[i];
      Point nrml = vec0 ^ vec1;
      nrml /= norm(nrml);

      double dp = nrml * normal_;
      if (std::abs(dp - 1.0) > TOL) {
        Errors::Message mesg("Polygon region is not exactly planar");
        Exceptions::amanzi_throw(mesg);
      }
    }

    /* Determine which direction to eliminate while projecting to one
       of the coordinate planes - to do this we have to find the
       direction in which the polygon is the smallest or in other
       words, the direction in which the normal to the polygon is the
       largest */
    int dmax = -1;
    double maxlen = -1;
    for (int i = 0; i < 3; i++) {
      double tmp = std::fabs(normal_[i]);
      if (tmp > maxlen) {
        maxlen = tmp;
        dmax = i;
      }
    }

    elim_dir_ = dmax;

  } else {
    Errors::Message mesg;
    mesg << "Cannot handle polygon regions with points of dimension " << space_dimension;
    Exceptions::amanzi_throw(mesg);
  }
}

// -------------------------------------------------------------
// RegionPolygon::inside
// -------------------------------------------------------------
bool
RegionPolygon::inside(const Point& p) const
{
#ifdef ENABLE_DBC
  if (p.dim() != points_[0].dim()) {
    Errors::Message mesg;
    mesg << "Mismatch in corner dimension of Polygon \"" << get_name() << "\" and query point.";
    Exceptions::amanzi_throw(mesg);
  }
#endif

  /* First check if the point is on the infinite line/plane */
  double d(0.0), res(0.0);

  for (int i = 0; i < p.dim(); ++i) {
    res += normal_[i] * p[i];
    d += normal_[i] * points_[0][i];
  }
  res -= d;

  if (std::abs(res) > TOL) return false;

  bool result(false);
  if (points_[0].dim() == 2) {
    // Now check if it lies in the line segment

    // vector from start of segment to point
    Point vec0 = p - points_[0];

    // segment vector
    Point vec1 = points_[1] - points_[0];

    // Normalize
    double slen = norm(vec1);
    vec1 /= slen;

    double dp = vec0 * vec1;

    // projection of vec0 along segment lies inside the segment
    if (dp >= 0 && dp <= slen) {
      // projected point along line segment
      Point p1 = points_[0] + dp * vec1;

      // vector between point and its projection
      Point dvec = p - p1;

      // distance between point and its projection
      double d_sqr = L22(dvec);

      // Is the distance 0? Point is inside segment
      if (d_sqr < TOL * TOL) result = true;
    }

  } else {
    /* Now check if the point is in the polygon */

    /* 
       Find the indices of coordinates on the projection plane
       
       if elim_dir_ is 0, then d0 = 1, d1 = 2 (YZ plane)
       if elim_dir_ is 1, then d0 = 2, d1 = 0 (XZ plane)
       if elim_dir_ is 2, then d0 = 0, d1 = 1 (XY plane)
    */

    double d0 = (elim_dir_ + 1) % 3;
    double d1 = (elim_dir_ + 2) % 3;

    /* Now apply the Jordan curve theorem to do the in/out test */
    /* odd number of crossings - point is inside                */

    double u, v;
    u = p[d0];
    v = p[d1];

    for (int i = 0; i != size(); ++i) {
      int iplus1 = (i + 1) % size();
      double u_i = points_[i][d0];
      double v_i = points_[i][d1];
      double u_iplus1 = points_[iplus1][d0];
      double v_iplus1 = points_[iplus1][d1];

      // don't compute - v_iplus1-v_i could be zero
      //      double slope = (u_iplus1-u_i)/(v_iplus1-v_i);

      if (((v_i > v && v_iplus1 <= v) || (v_iplus1 > v && v_i <= v)) &&
          (u <= (u_i + (v - v_i) * (u_iplus1 - u_i) / (v_iplus1 - v_i))))
        result = !result;
    }

    // The above check is not guaranteed to give an +ve result if the point is on
    // the boundary. So do an additional check

    if (!result) {
      for (int i = 0; i < size(); i++) {
        int iplus1 = (i + 1) % size();

        Point p_i(2);
        p_i.set(points_[i][d0], points_[i][d1]);

        Point p_iplus1(2);
        p_iplus1.set(points_[iplus1][d0], points_[iplus1][d1]);

        // vector from first point of segment to query point
        Point vec0(2);
        vec0.set(u - p_i[0], v - p_i[1]);

        // line segment vector
        Point vec1(2);
        vec1 = p_iplus1 - p_i;

        // unit vector along segment
        double slen = norm(vec1);
        vec1 /= slen;

        double dp = vec0 * vec1;

        // projection of vec0 along segment lies outside the segment
        if (dp < 0 || dp > slen) continue;

        // projected point along line segment
        Point p1(2);
        p1 = p_i + dp * vec1;

        // vector between projected point and query point
        Point dvec(2);
        dvec.set(u - p1[0], v - p1[1]);

        // distance between point and its projection
        double d_sqr = L22(dvec);

        // Is the distance 0? Point is inside segment
        if (d_sqr < TOL * TOL) {
          result = true;
          break;
        }
      }
    }
  }

  return result;
}


} // namespace AmanziGeometry
} // namespace Amanzi
