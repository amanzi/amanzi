/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "Geometry.hh"

using namespace std;

namespace Amanzi {
namespace AmanziGeometry {

// Checks if point is inside polyhedron
//
// ccoords  - vertices of the polyhedron (in no particular order)
// nf       - number of faces of polyhedron
// nfnodes  - number of nodes for each face
// fcoords  - linear array of face coordinates in in ccw manner
//            assuming normal of face is pointing out (
//
// So if the polyhedron has 5 faces with 5,3,3,3 and 3 nodes each
// then entries 1-5 of fcoords describes face 1, entries 6-9
// describes face 2 and so on
//
// Assuming that the polyhedron's faces can be broken into
// triangular subfaces, this routine checks that the test point
// forms a positive volume with each triangular subface
bool
point_in_polyhed(const Point testpnt,
                 const std::vector<Point>& ccoords,
                 const std::size_t nf,
                 const std::vector<std::size_t>& nfnodes,
                 const std::vector<Point>& fcoords)
{
  int np = ccoords.size();
  if (np < 4) {
    cout << "Not a polyhedron" << std::endl;
    return false;
  }

  int offset = 0;
  for (int i = 0; i < nf; i++) {
    Point v1(3), v2(3), v3(3);
    double tvolume;

    if (nfnodes[i] == 3) {
      v1 = fcoords[offset] - testpnt;
      v2 = fcoords[offset + 1] - testpnt;
      v3 = fcoords[offset + 2] - testpnt;
      tvolume = (v1 ^ v2) * v3;

      if (tvolume < 0.0) return false;
    } else {
      // geometric center of all face nodes
      Point fcenter(0.0, 0.0, 0.0);
      for (int j = 0; j < nfnodes[i]; j++) fcenter += fcoords[offset + j];
      fcenter /= nfnodes[i];

      for (int j = 0; j < nfnodes[i]; j++) { // for each edge of face
        // form tet from edge of face, face center and test point
        Point tcentroid(3);
        int k, kp1;

        k = offset + j;
        kp1 = offset + (j + 1) % nfnodes[i];

        v1 = fcoords[k] - testpnt;
        v2 = fcoords[kp1] - testpnt;
        v3 = fcenter - testpnt;
        tvolume = (v1 ^ v2) * v3;

        if (tvolume < 0.0) return false;
      } // for each edge of face

      offset += nfnodes[i];
    }
  } // for each face

  return true;
} // point_in_polyhed


// Check if point is in polygon by Jordan's crossing algorithm
bool
point_in_polygon(const Point testpnt, const std::vector<Point>& coords)
{
  int i, ip1, c;

  /* Basic test - will work for strictly interior and exterior points */
  int np = coords.size();

  double x = testpnt.x();
  double y = testpnt.y();

  for (i = 0, c = 0; i < np; i++) {
    //    std::cout<<"coords "<<coords[i]<<"\n";
    ip1 = (i + 1) % np;
    if (((coords[i].y() > y && coords[ip1].y() <= y) ||
         (coords[ip1].y() > y && coords[i].y() <= y)) &&
        (x <= (coords[i].x() + (y - coords[i].y()) * (coords[ip1].x() - coords[i].x()) /
                                 (coords[ip1].y() - coords[i].y()))))
      c = !c;
  }

  /* If we don't need consistent classification of points on the
     boundary, we can quit here.  If the point is classified as
     inside, it is definitely inside or on the boundary - no way it
     can be outside and be classified as inside */

  return (c == 1);
}

} // namespace AmanziGeometry
} // namespace Amanzi
