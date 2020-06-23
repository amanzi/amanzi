/*
  Geometry

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Lipnikov Konstantin (lipnikov@lanl.gov)
           Rao Garimella (rao@lanl.gov)

  A box region defined by two corner points and normals to its side.
*/

#include <vector>

#include "dbc.hh"
#include "errors.hh"

#include "Point.hh"
#include "RegionBoxVolumeFractions.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------------
// Constructor
// -------------------------------------------------------------------
RegionBoxVolumeFractions::RegionBoxVolumeFractions(
    const std::string& name, const int id,
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


// -------------------------------------------------------------------
// Implementation of a virtual member function.
// -------------------------------------------------------------------
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


// -------------------------------------------------------------------
// Implementation of a virtual member function.
// We have to analyze 
// -------------------------------------------------------------------
double RegionBoxVolumeFractions::intersect(
    const std::vector<Point>& polytope,
    const std::vector<std::vector<int> >& faces) const
{
  double volume(0.0);
  int sdim = polytope[0].dim();

  if ((sdim == 2 && degeneracy_ < 0) || (sdim == 3 && degeneracy_ >= 0)) {
    std::vector<Point> box, result_xy;

    box.push_back(Point(0.0, 0.0));
    box.push_back(Point(1.0, 0.0)); 
    box.push_back(Point(1.0, 1.0));
    box.push_back(Point(0.0, 1.0));

    // project manifold on a plane
    Point p3d(3), p2d(2);
    std::vector<Point> nodes;

    if (degeneracy_ >= 0) {
      double eps(1.0e-6);

      for (int n = 0; n < polytope.size(); ++n) {
        p3d = N_ * (polytope[n] - p0_);
        if (std::abs(p3d[degeneracy_]) > eps) return volume;

        int i(0);
        for (int k = 0; k < sdim; ++k) p2d[i++] = p3d[k];
        nodes.push_back(p2d);
      }
      IntersectConvexPolygons(nodes, box, result_xy);
    } else {
      for (int n = 0; n < polytope.size(); ++n) {
        p2d = N_ * (polytope[n] - p0_);
        nodes.push_back(p2d);
      }
    }

    IntersectConvexPolygons(nodes, box, result_xy);

    int nnodes = result_xy.size(); 
    if (nnodes > 0) {
      for (int i = 0; i < nnodes; ++i) {
        int j = (i + 1) % nnodes;

        const Point& p1 = result_xy[i];
        const Point& p2 = result_xy[j];
        volume += (p1[0] + p2[0]) * (p2[1] - p1[1]) / 2;
      }
      volume *= jacobian_;
    }
  }

  else if (sdim == 3 && degeneracy_ < 0) {
    std::vector<std::vector<int> > result_faces;
    std::vector<Point> result_xyz;
    std::vector<std::pair<Point, Point> > box;

    box.push_back(std::make_pair(Point(0.0, 0.0, 0.0), Point(-1.0, 0.0, 0.0)));
    box.push_back(std::make_pair(Point(0.0, 0.0, 0.0), Point(0.0, -1.0, 0.0)));
    box.push_back(std::make_pair(Point(0.0, 0.0, 0.0), Point(0.0, 0.0, -1.0)));

    box.push_back(std::make_pair(Point(1.0, 1.0, 1.0), Point(1.0, 0.0, 0.0)));
    box.push_back(std::make_pair(Point(1.0, 1.0, 1.0), Point(0.0, 1.0, 0.0)));
    box.push_back(std::make_pair(Point(1.0, 1.0, 1.0), Point(0.0, 0.0, 1.0)));

    std::vector<Point> nodes;
    Point p3d(3);

    for (int n = 0; n < polytope.size(); ++n) {
      p3d = N_ * (polytope[n] - p0_);
      nodes.push_back(p3d);
    }
    IntersectConvexPolyhedra(nodes, faces, box, result_xyz, result_faces);

    int nfaces = result_faces.size(); 
    if (nfaces > 3) {
      for (int i = 0; i < nfaces; ++i) {
        int nnodes = result_faces[i].size();
        for (int k = 0; k < nnodes - 2; ++k) {
          const Point& p0 = result_xyz[result_faces[i][k]];
          const Point& p1 = result_xyz[result_faces[i][k + 1]];
          const Point& p2 = result_xyz[result_faces[i][k + 2]];
          volume += ((p1 - p0)^(p2 - p0)) * p0;
        }
      }
      volume *= jacobian_ / 6;
    }
  } else {
    AMANZI_ASSERT(0);
  }

  return volume;
}


// -------------------------------------------------------------------
// Non-member function.
// Intersection of two counter clockwise oriented polygons given
// by vertices xy1 and xy2.
// -------------------------------------------------------------------
void IntersectConvexPolygons(const std::vector<Point>& xy1,
                             const std::vector<Point>& xy2,
                             std::vector<Point>& xy3)
{
  std::list<std::pair<double, Point> > result;
  std::list<std::pair<double, Point> >::iterator it, it_next, it2;

  // populate list with the second polygon
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


// -------------------------------------------------------------------
// Non-member function.
// Intersection of two convex polyhedra P1 and P2. Polyhedron P1 is 
// defined by vertices xyz1 and faces ordered counter clockwise with
// respect to their exterior normals. Polyhedron P2 is defined by a
// set of half-spaces (point and exterior normal). The result is
// polyhedron P3 ordered similar to P1.
// -------------------------------------------------------------------
// #define VERBOSE
void IntersectConvexPolyhedra(const std::vector<Point>& xyz1,
                              const std::vector<std::vector<int> >& faces1,
                              const std::vector<std::pair<Point, Point> >& xyz2,
                              std::vector<Point>& xyz3,
                              std::vector<std::vector<int> >& faces3)
{
  // initialize dynamic polyhedron structure with the first polyhedron
  int nfaces1 = faces1.size();
  std::vector<std::pair<double, Point> > result_xyz;
  std::list<ClippedFace> result(nfaces1);

  int nxyz = xyz1.size();
  for (int i = 0; i < nxyz; ++i) {
    result_xyz.push_back(std::make_pair(0.0, xyz1[i]));
  }

  int k(0);
  std::list<ClippedFace>::iterator itf;

  for (itf = result.begin(); itf != result.end(); ++itf, ++k) {
    for (int n = 0; n < faces1[k].size(); ++n) {
      itf->nodes.push_back(faces1[k][n]);
    }
  }

  // clip polyhedron using the second polyhedron.
  int nfaces2 = xyz2.size();
  double d1, d2, tmp, eps(1e-7);
  Point v1(xyz1[0]);
  std::list<int>::iterator itv, itv_prev, itv_next, itv2;

  for (int n = 0; n < nfaces2; ++n) {
    const Point& p = xyz2[n].first;
    const Point& normal = xyz2[n].second;
#ifdef VERBOSE
std::cout << "plane=" << p << " normal=" << normal << std::endl;
#endif

    // location of nodes relative to the n-th plane
    for (int i = 0; i < nxyz; ++i) {
      if (result_xyz[i].first <= 0.0) {
        tmp = normal * (result_xyz[i].second - p);
        if (std::fabs(tmp) < eps) tmp = 0.0;
        result_xyz[i].first = tmp / norm(normal);
      }
    }

    // clip each face of the clipped polyhedron "result"
    if (result.size() == 0) continue;

    for (itf = result.begin(); itf != result.end(); ++itf) {
      itf->new_edge = std::make_pair(-1, -1);
      itf->edge_flag = 0;
    }

    for (itf = result.begin(); itf != result.end(); ++itf) {
      for (itv = itf->nodes.begin(); itv != itf->nodes.end(); ++itv) {
        itv_next = itv;
        if (++itv_next == itf->nodes.end()) itv_next = itf->nodes.begin();

        std::pair<double, Point>& p1 = result_xyz[*itv];
        std::pair<double, Point>& p2 = result_xyz[*itv_next];

        d1 = p1.first;
        d2 = p2.first;
        // add vertex if intersection was found; otherwise, remove vertex.
        if (d1 * d2 < 0.0) {  
          tmp = d2 / (d2 - d1);
          v1 = tmp * p1.second + (1.0 - tmp) * p2.second; 

          int idx(-1);
          for (int i = 0; i < nxyz; ++i) {
            if (norm(v1 - result_xyz[i].second) < 1e-8) {
              idx = i;
              break;
            }
          }

          if (idx >= 0) {
            itf->nodes.insert(itv_next, idx);
          } else {
            result_xyz.push_back(std::make_pair(0.0, v1));
            itf->nodes.insert(itv_next, nxyz);
#ifdef VERBOSE
std::cout << "  add " << v1 << std::endl;
#endif
            nxyz++;
          }
        }
        // this edge may be nedeed for building new face 
        else if (d1 == 0.0 && d2 == 0.0) {
          itf->new_edge = std::make_pair(*itv, *itv_next);
          itf->edge_flag = 0;
        }
      }

      // removing cut-out edges for each face separately
      itv = itf->nodes.begin();
      while (itv != itf->nodes.end()) {
        if (result_xyz[*itv].first > 0.0) {
          itv = itf->nodes.erase(itv);
#ifdef VERBOSE
if (itf->nodes.size() == 2) {
std::cout << "  removing face: ";
for (std::list<int>::iterator itt = itf->nodes.begin(); itt != itf->nodes.end(); ++itt) std::cout << *itt << " ";
std::cout << std::endl;
}
#endif
          if (itf->nodes.size() == 2) {
            itf = result.erase(itf);
            itf--;
            break;
          }

          itv_next = itv;
          if (itv_next == itf->nodes.end()) itv_next = itf->nodes.begin();

          itv_prev = itv;
          if (itv_prev == itf->nodes.begin()) itv_prev = itf->nodes.end();
          else itv_prev--;

          itf->new_edge = std::make_pair(*itv_prev, *itv_next);
          itf->edge_flag = 1;
#ifdef VERBOSE
std::cout << "  new edge:" << *itv_prev << " " << *itv_next << "   p1=" 
          << result_xyz[*itv_prev].second << "  p2=" << result_xyz[*itv_next].second 
          << " d=" << result_xyz[*itv_prev].first << " " << result_xyz[*itv_next].first << std::endl;
#endif
        } else {
          itv++;
        }
      }
    } 
    if (result.size() < 4) continue;

    // forming a new face
    // -- we have enough new edges
    int nedges(0), edge_flag(0);
    for (itf = result.begin(); itf != result.end(); ++itf) {
      if (itf->new_edge.second >= 0) nedges++;
      edge_flag = std::max(edge_flag, itf->edge_flag);
    }

    // -- starting point for the new face
    if (nedges > 2 && edge_flag == 1) {
#ifdef VERBOSE
std::cout << "  enough edges to build new face: " << nedges << std::endl;
#endif
      int n1, n2, n3, n4(-1);
      ClippedFace new_face;

      for (itf = result.begin(); itf != result.end(); ++itf) {
        n1 = itf->new_edge.second;
        n2 = itf->new_edge.first;
        if (n1 >= 0) {
          new_face.nodes.push_back(n1);
          break;
        }
      }

      // -- the remaining points of the new face
      while(n4 != n1) {
        for (itf = result.begin(); itf != result.end(); ++itf) {
          n3 = itf->new_edge.second;
          n4 = itf->new_edge.first;
          if (n2 == n3) {
            new_face.nodes.push_back(n3);
            n2 = n4;
            break;
          }
        }
      }

      if (new_face.nodes.size() > 2) {
#ifdef VERBOSE
std::cout << "  adding new face nodes: ";
for (itv = new_face.nodes.begin(); itv != new_face.nodes.end(); ++itv) std::cout << *itv << " ";
std::cout << std::endl;
#endif
        result.push_back(new_face);
      }
    }
  }

  // output of the result
  xyz3.clear();
  faces3.clear();

  if (result.size() > 3) { 
    int nxyz3(result_xyz.size());
    std::vector<int> map(nxyz3, -1);

    faces3.resize(result.size());

    k = 0;
    for (itf = result.begin(); itf != result.end(); ++itf, ++k) {
      for (itv = itf->nodes.begin(); itv != itf->nodes.end(); ++itv) {
        faces3[k].push_back(*itv);
        map[*itv] = 0;
      }
    }

    int m(0);
    for (int i = 0; i < nxyz3; ++i) {
      if (map[i] == 0) map[i] = m++;
    }

    // -- count true vertices
    for (int i = 0; i < nxyz3; ++i) {
      if (map[i] >= 0) xyz3.push_back(result_xyz[i].second);
    }

#ifdef VERBOSE
std::cout << "updating map" << std::endl;
#endif
    // -- update face-to-nodes maps
    for (int i = 0; i < result.size(); ++i) {
      int nnodes = faces3[i].size();
      for (int n = 0; n < nnodes; ++n) {
        faces3[i][n] = map[faces3[i][n]];
      }
    }
  }
}

}  // namespace AmanziGeometry
}  // namespace Amanzi
