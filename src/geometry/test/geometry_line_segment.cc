/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <UnitTest++.h>

#include <cmath>
#include <vector>

#include "Point.hh"
#include "RegionLineSegment.hh"

using namespace Amanzi::AmanziGeometry;

const double tol = 1e-10;

namespace {

// convex pentagon with ordered (counter-clockwise) vertices
std::vector<Point>
Pentagon()
{
  return { Point(0.0, 0.0), Point(4.0, 0.0), Point(5.0, 2.0), Point(2.5, 4.0), Point(-1.0, 2.0) };
}


// general parallelepiped P(i,j,k) = i * (2,0,0) + j * (1,2,0) + k * (0.5,0.5,2),
// 0 <= i,j,k <= 1; vertex number is i + 2 * j + 4 * k
std::vector<Point>
Parallelepiped(std::vector<std::vector<int>>& faces)
{
  faces = { { 0, 2, 3, 1 }, { 4, 5, 7, 6 }, { 0, 1, 5, 4 }, { 2, 6, 7, 3 }, { 0, 4, 6, 2 }, { 1, 3, 7, 5 } };
  return { Point(0.0, 0.0, 0.0), Point(2.0, 0.0, 0.0), Point(1.0, 2.0, 0.0), Point(3.0, 2.0, 0.0),
           Point(0.5, 0.5, 2.0), Point(2.5, 0.5, 2.0), Point(1.5, 2.5, 2.0), Point(3.5, 2.5, 2.0) };
}

} // namespace



// ---------------------------------------------------------------
// 2D tests
// ---------------------------------------------------------------
TEST(LINE_SEGMENT_2D_REGULAR)
{
  auto polygon = Pentagon();
  std::vector<std::vector<int>> faces; // not used in 2D

  { // crosses two edges: y = 1 between x = -0.5 and x = 4.5
    RegionLineSegment segment("segment", 1, Point(-3.0, 1.0), Point(7.0, 1.0));
    CHECK_CLOSE(5.0, segment.intersect(polygon, faces), tol);
    Point mid(2);
    segment.ComputeInterLinePoints(polygon, faces, mid);
    CHECK_CLOSE(0.0, norm(mid - Point(2.0, 1.0)), tol);
  }
  { // completely inside
    RegionLineSegment segment("segment", 1, Point(1.0, 1.0), Point(3.0, 2.0));
    CHECK_CLOSE(std::sqrt(5.0), segment.intersect(polygon, faces), tol);
    Point mid(2);
    segment.ComputeInterLinePoints(polygon, faces, mid);
    CHECK_CLOSE(0.0, norm(mid - Point(2.0, 1.5)), tol);
  }
  { // one end inside, exits through the edge (2.5,4)-(-1,2) at y = 26/7
    RegionLineSegment segment("segment", 1, Point(2.0, 1.0), Point(2.0, 6.0));
    CHECK_CLOSE(19.0 / 7.0, segment.intersect(polygon, faces), tol);
    Point mid(2);
    segment.ComputeInterLinePoints(polygon, faces, mid);
    CHECK_CLOSE(0.0, norm(mid - Point(2.0, 0.5 * (1.0 + 26.0 / 7.0))), tol);
  }
  { // the same, but the first end is outside
    RegionLineSegment segment("segment", 1, Point(2.0, 6.0), Point(2.0, 1.0));
    CHECK_CLOSE(19.0 / 7.0, segment.intersect(polygon, faces), tol);
  }
}


TEST(LINE_SEGMENT_2D_NO_INTERSECTION)
{
  auto polygon = Pentagon();
  std::vector<std::vector<int>> faces; // not used in 2D

  { // line misses the polygon
    RegionLineSegment segment("segment", 1, Point(6.0, 0.0), Point(6.0, 5.0));
    CHECK_CLOSE(0.0, segment.intersect(polygon, faces), tol);
  }
  { // line crosses the polygon but the segment is on the left
    RegionLineSegment segment("segment", 1, Point(-4.0, 1.0), Point(-3.0, 1.0));
    CHECK_CLOSE(0.0, segment.intersect(polygon, faces), tol);
  }
  { // line crosses the polygon but the segment is on the right
    RegionLineSegment segment("segment", 1, Point(6.0, 1.0), Point(8.0, 1.0));
    CHECK_CLOSE(0.0, segment.intersect(polygon, faces), tol);
  }
}


TEST(LINE_SEGMENT_2D_VERTICES_AND_EDGES)
{
  auto polygon = Pentagon();
  std::vector<std::vector<int>> faces; // not used in 2D

  { // goes through vertex (4,0) and exits via the edge (5,2)-(2.5,4) at y = 2.8
    RegionLineSegment segment("segment", 1, Point(4.0, -1.0), Point(4.0, 5.0));
    CHECK_CLOSE(2.8, segment.intersect(polygon, faces), tol);
    Point mid(2);
    segment.ComputeInterLinePoints(polygon, faces, mid);
    CHECK_CLOSE(0.0, norm(mid - Point(4.0, 1.4)), tol);
  }
  { // goes through two vertices, (0,0) and (5,2)
    RegionLineSegment segment("segment", 1, Point(-5.0, -2.0), Point(10.0, 4.0));
    CHECK_CLOSE(std::sqrt(29.0), segment.intersect(polygon, faces), tol);
    Point mid(2);
    segment.ComputeInterLinePoints(polygon, faces, mid);
    CHECK_CLOSE(0.0, norm(mid - Point(2.5, 1.0)), tol);
  }
  { // lies on the edge (0,0)-(4,0)
    RegionLineSegment segment("segment", 1, Point(-1.0, 0.0), Point(5.0, 0.0));
    CHECK_CLOSE(4.0, segment.intersect(polygon, faces), tol);
  }
  { // ends at a vertex from inside
    RegionLineSegment segment("segment", 1, Point(2.0, 1.0), Point(5.0, 2.0));
    CHECK_CLOSE(std::sqrt(10.0), segment.intersect(polygon, faces), tol);
  }
}


TEST(LINE_SEGMENT_2D_CLOCKWISE_POLYGON)
{
  // the same pentagon with clockwise vertices; the vertex order must not matter
  auto ccw = Pentagon();
  std::vector<Point> polygon(ccw.rbegin(), ccw.rend());
  std::vector<std::vector<int>> faces; // not used in 2D

  { // crosses two edges
    RegionLineSegment segment("segment", 1, Point(-3.0, 1.0), Point(7.0, 1.0));
    CHECK_CLOSE(5.0, segment.intersect(polygon, faces), tol);
  }
  { // goes through a vertex
    RegionLineSegment segment("segment", 1, Point(4.0, -1.0), Point(4.0, 5.0));
    CHECK_CLOSE(2.8, segment.intersect(polygon, faces), tol);
  }
  { // line misses the polygon
    RegionLineSegment segment("segment", 1, Point(6.0, 0.0), Point(6.0, 5.0));
    CHECK_CLOSE(0.0, segment.intersect(polygon, faces), tol);
  }
}


// ---------------------------------------------------------------
// 3D tests
// ---------------------------------------------------------------
TEST(LINE_SEGMENT_3D_REGULAR)
{
  std::vector<std::vector<int>> faces;
  auto polytope = Parallelepiped(faces);

  { // crosses two opposite faces through their centers
    RegionLineSegment segment("segment", 1, Point(-0.25, 1.25, 1.0), Point(3.75, 1.25, 1.0));
    CHECK_CLOSE(2.0, segment.intersect(polytope, faces), tol);
    Point mid(3);
    segment.ComputeInterLinePoints(polytope, faces, mid);
    CHECK_CLOSE(0.0, norm(mid - Point(1.75, 1.25, 1.0)), tol);
  }
  { // completely inside
    RegionLineSegment segment("segment", 1, Point(0.9, 0.8, 0.8), Point(2.25, 1.45, 1.0));
    CHECK_CLOSE(std::sqrt(2.285), segment.intersect(polytope, faces), tol);
    Point mid(3);
    segment.ComputeInterLinePoints(polytope, faces, mid);
    CHECK_CLOSE(0.0, norm(mid - Point(1.575, 1.125, 0.9)), tol);
  }
  { // one end inside, exits through the top face at (2,1.5,2)
    RegionLineSegment segment("segment", 1, Point(1.75, 1.25, 1.0), Point(2.25, 1.75, 3.0));
    CHECK_CLOSE(std::sqrt(1.125), segment.intersect(polytope, faces), tol);
    Point mid(3);
    segment.ComputeInterLinePoints(polytope, faces, mid);
    CHECK_CLOSE(0.0, norm(mid - Point(1.875, 1.375, 1.5)), tol);
  }
  { // the same, but the first end is outside
    RegionLineSegment segment("segment", 1, Point(2.25, 1.75, 3.0), Point(1.75, 1.25, 1.0));
    CHECK_CLOSE(std::sqrt(1.125), segment.intersect(polytope, faces), tol);
  }
}


TEST(LINE_SEGMENT_3D_NO_INTERSECTION)
{
  std::vector<std::vector<int>> faces;
  auto polytope = Parallelepiped(faces);

  { // line misses the parallelepiped
    RegionLineSegment segment("segment", 1, Point(-2.75, -1.75, 1.0), Point(0.25, 4.25, 1.0));
    CHECK_CLOSE(0.0, segment.intersect(polytope, faces), tol);
  }
  { // line crosses the parallelepiped but the segment is on the left
    RegionLineSegment segment("segment", 1, Point(-1.25, 1.25, 1.0), Point(-0.25, 1.25, 1.0));
    CHECK_CLOSE(0.0, segment.intersect(polytope, faces), tol);
  }
  { // line crosses the parallelepiped but the segment is on the right
    RegionLineSegment segment("segment", 1, Point(3.75, 1.25, 1.0), Point(4.75, 1.25, 1.0));
    CHECK_CLOSE(0.0, segment.intersect(polytope, faces), tol);
  }
  { // segment is parallel to the face 4-5-7-6 (z = 2) and lies above the parallelepiped
    RegionLineSegment segment("segment", 1, Point(0.0, 0.0, 3.0), Point(4.0, 3.0, 3.0));
    CHECK_CLOSE(0.0, segment.intersect(polytope, faces), tol);
  }
  { // segment is parallel to the face 0-2-3-1 (z = 0) and lies below the parallelepiped
    RegionLineSegment segment("segment", 1, Point(0.0, 1.0, -1.0), Point(3.0, 1.0, -1.0));
    CHECK_CLOSE(0.0, segment.intersect(polytope, faces), tol);
  }
  { // segment is parallel to the faces 0-1-5-4 and 2-6-7-3 (normal is (0,-4,1)) and
    // lies behind the face 2-6-7-3
    RegionLineSegment segment("segment", 1, Point(0.0, 4.0, 4.0), Point(4.0, 4.0, 4.0));
    CHECK_CLOSE(0.0, segment.intersect(polytope, faces), tol);
  }
}


TEST(LINE_SEGMENT_3D_VERTICES_AND_EDGES)
{
  std::vector<std::vector<int>> faces;
  auto polytope = Parallelepiped(faces);

  { // main diagonal: goes through vertices 0 and 7
    RegionLineSegment segment("segment", 1, Point(-1.75, -1.25, -1.0), Point(5.25, 3.75, 3.0));
    CHECK_CLOSE(std::sqrt(22.5), segment.intersect(polytope, faces), tol);
    Point mid(3);
    segment.ComputeInterLinePoints(polytope, faces, mid);
    CHECK_CLOSE(0.0, norm(mid - Point(1.75, 1.25, 1.0)), tol);
  }
  { // crosses the edge 0-1 at (1,0,0) and ends inside
    RegionLineSegment segment("segment", 1, Point(0.25, -1.25, -1.0), Point(1.75, 1.25, 1.0));
    CHECK_CLOSE(std::sqrt(3.125), segment.intersect(polytope, faces), tol);
  }
  { // enters through the edge 0-1 at (1,0,0) and exits through the edge 6-7 at (2.5,2.5,2)
    RegionLineSegment segment("segment", 1, Point(0.25, -1.25, -1.0), Point(3.25, 3.75, 3.0));
    CHECK_CLOSE(std::sqrt(12.5), segment.intersect(polytope, faces), tol);
    Point mid(3);
    segment.ComputeInterLinePoints(polytope, faces, mid);
    CHECK_CLOSE(0.0, norm(mid - Point(1.75, 1.25, 1.0)), tol);
  }
  { // starts at the vertex 0 and goes inside
    RegionLineSegment segment("segment", 1, Point(0.0, 0.0, 0.0), Point(1.75, 1.25, 1.0));
    CHECK_CLOSE(std::sqrt(5.625), segment.intersect(polytope, faces), tol);
  }
}
