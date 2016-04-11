/*
  A rectangular region in space, defined by two corner points and
  normals to its side.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Lipnikov Konstantin (lipnikov@lanl.gov)
           Rao Garimella (rao@lanl.gov)
*/

#include <vector>

#include "dbc.hh"
#include "errors.hh"

#include "Point.hh"
#include "RegionBoxVolumeFractions.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
RegionBoxVolumeFractions::RegionBoxVolumeFractions(
    const std::string& name, const Set_ID id,
    const Point& p0, const Point& p1,
    const std::vector<Point>& normals,
    const LifeCycleType lifecycle)
  : Region(name, id, true, BOX_VOF, p0.dim(), p0.dim(), lifecycle),
    p0_(p0),
    p1_(p1),
    normals_(normals),
    degeneracy_(-1)
{
  Errors::Message msg;
  if (p0_.dim() != p1_.dim()) {
    msg << "Mismatch in dimensions of corner points of RegionBoxVolumeFractions \""
        << Region::name() << "\"";
    Exceptions::amanzi_throw(msg);

    for (int n = 0; n < normals.size(); ++n) {
      if (p0_.dim() != normals_[n].dim()) {
        msg << "Mismatch in dimensions of points and normals of RegionBoxVolumeFractions \""
            << Region::name() << "\"";
        Exceptions::amanzi_throw(msg);
      }
    }
  }

  // Calculate the transformation tensor
  int dim = p0.dim();
  N_.set(dim);

  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      N_(j, i) = normals_[i][j];
    }
  }
  N_.Inverse();

  // Check if this is a reduced dimensionality box (e.g. even though
  // it is in 3D space it is a 2D box). Normalize the tensor to unit
  // reference cube.
  Point v1(dim);
  v1 = N_ * (p1_ - p0_);

  jacobian_ = 1.0;
  int mdim(dim);
  for (int i = 0; i != dim; ++i) {
    double len = v1[i];
    if (std::abs(len) < 1e-12 * norm(v1)) {
      mdim--;
      degeneracy_ = i;
    } else {
      jacobian_ *= std::abs(len);
      for (int j = 0; j < dim; ++j) N_(i, j) /= len;
    }
  } 
  
  if (mdim < dim) set_manifold_dimension(mdim);

  if (mdim < dim - 1) {
    msg << "The box degenerate in more than one direction is not supported.\n";
    Exceptions::amanzi_throw(msg);
  }
}


// -------------------------------------------------------------
// Implementation of a virtual member function.
// -------------------------------------------------------------
bool RegionBoxVolumeFractions::inside(const Point& p) const
{
#ifdef ENABLE_DBC
  if (p.dim() != p0_.dim()) {
    Errors::Message msg;
    msg << "Mismatch in corner dimension of RegionBox \""
        << name() << "\" and query point.";
    Exceptions::amanzi_throw(msg);
  }
#endif

  Point phat(N_ * (p - p0_));
 
  bool result(true);
  for (int i = 0; i != p.dim(); ++i) {
    result &= (phat[i] > -TOL && phat[i] < 1.0 + TOL);
  }
  return result;
}


// -------------------------------------------------------------
// Implementation of a virtual member function.
// We have to analyze 
// -------------------------------------------------------------
double RegionBoxVolumeFractions::intersect(const std::vector<Point>& polytope) const
{
  double volume(-1.0);
  int mdim, sdim;
  std::vector<Point> box, result;

  mdim = manifold_dimension();
  sdim = polytope[0].dim();
   
  if ((sdim == 2 && degeneracy_ < 0) || (sdim == 3 && degeneracy_ >= 0)) {
    box.push_back(Point(0.0, 0.0));
    box.push_back(Point(1.0, 0.0)); 
    box.push_back(Point(1.0, 1.0));
    box.push_back(Point(0.0, 1.0));

    // project manifold on a plane
    if (degeneracy_ >= 0) {
      double eps(1.0e-6);
      Point p3d(3), p2d(2);
      std::vector<Point> nodes;

      for (int n = 0; n < polytope.size(); ++n) {
        p3d = N_ * (polytope[n] - p0_);
        if (std::abs(p3d[degeneracy_]) > eps) return volume;

        int i(0);
        for (int k = 0; k < sdim; ++k) p2d[i++] = p3d[k];
        nodes.push_back(p2d);
      }
      IntersectConvexPolygons(nodes, box, result);
    } else {
      IntersectConvexPolygons(polytope, box, result);
    }

    int nnodes = result.size(); 
    if (nnodes > 0) {
      volume = 0.0;
      for (int i = 0; i < nnodes; ++i) {
        int j = (i + 1) % nnodes;

        const Point& p1 = result[i];
        const Point& p2 = result[j];
        volume += (p1[0] + p2[0]) * (p2[1] - p1[1]) / 2;
      }

      volume *= jacobian_;
    }
  }

  return volume;
}


// -------------------------------------------------------------
// Non-member function.
// Intersection of two counter clockwise oriented polygons given
// by vertices xy1 and xy2.
// -------------------------------------------------------------
void IntersectConvexPolygons(const std::vector<Point>& xy1,
                             const std::vector<Point>& xy2,
                             std::vector<Point>& xy3)
{
  std::list<std::pair<double, Point> > result;
  std::list<std::pair<double, Point> >::iterator it, it_next, it2;

  // populate list with the second polygon
  int n2 = xy2.size();
  for (int i = 0; i < xy2.size(); ++i) {
    result.push_back(std::make_pair(0.0, xy2[i]));
  }

  // intersect each edge of the first polygon with the result
  int n1 = xy1.size();
  double eps(1e-6);
  Point v1(xy1[0]);

  for (int i = 0; i < n1; ++i) {
    if (result.size() <= 2) break;

    int j = (i + 1) % n1;
    Point edge(xy1[j] - xy1[i]);
    Point normal(edge[1], -edge[0]);
    normal /= norm(normal);

    // Calculate distance of polygon vertices to the plane defined by 
    // the point xy1[i] and exterior normal normal.
    for (it = result.begin(); it != result.end(); ++it) {
      double tmp = normal * (it->second - xy1[i]);
      if (std::fabs(tmp) < eps) tmp = 0.0;
      it->first = tmp;
    }

    double d1, d2, tmp;
    for (it = result.begin(); it != result.end(); ++it) {
      it_next = it;
      if (++it_next == result.end()) it_next = result.begin();

      d1 = it->first;
      d2 = it_next->first;
      // add vertex if intersection was found; otherwise, remove vertex.
      if (d1 * d2 < 0.0) {  
        tmp = d2 / (d2 - d1);
        v1 = tmp * it->second + (1.0 - tmp) * it_next->second; 
        result.insert(it_next, std::make_pair(0.0, v1));
      }
    } 

    // removing cut-out edges
    it = result.begin();
    while (it != result.end()) {
      if (it->first > 0.0) {
        it = result.erase(it);
        if (result.size() == 2) break;
      } else {
        it++;
      }
    }
  }

  xy3.clear();
  if (result.size() > 2) { 
    for (it = result.begin(); it != result.end(); ++it) {
      xy3.push_back(it->second);
    }
  }
}

}  // namespace AmanziGeometry
}  // namespace Amanzi
