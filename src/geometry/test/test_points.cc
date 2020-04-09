/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#include <UnitTest++.h>

#include <iostream>


#include "../Point.hh"

#include "mpi.h"


TEST(Point)
{
  double x = 0.4, y = 0.6, z = 0.8;
  double xyz[3] = { 1, 2, 3.5 };

  // Create different types of points

  Amanzi::AmanziGeometry::Point p0(2);
  p0.set(x, y);
  CHECK_EQUAL(x, p0.x());
  CHECK_EQUAL(y, p0.y());

  Amanzi::AmanziGeometry::Point p1(2);
  p1.set(2, xyz);
  CHECK_EQUAL(xyz[0], p1.x());
  CHECK_EQUAL(xyz[1], p1.y());


  x = x + 0.1;
  y = y + 0.1;
  Amanzi::AmanziGeometry::Point p2(x, y);
  CHECK_EQUAL(x, p2.x());
  CHECK_EQUAL(y, p2.y());

  Amanzi::AmanziGeometry::Point p3(p2);
  CHECK_EQUAL(x, p3.x());
  CHECK_EQUAL(y, p3.y());

  Amanzi::AmanziGeometry::Point p4 = p2;
  CHECK_EQUAL(x, p4.x());
  CHECK_EQUAL(y, p4.y());

  Amanzi::AmanziGeometry::Point p5(3);
  xyz[0] = x;
  xyz[1] = y;
  xyz[2] = z;
  p5.set(3, xyz);
  CHECK_EQUAL(x, p5.x());
  CHECK_EQUAL(y, p5.y());
  CHECK_EQUAL(z, p5.z());


  // Create a bunch of 3D points and delete some of them to trigger
  // different patterns of deletion in the underlying coordinate
  // buffer

  int num3 = 25;
  Amanzi::AmanziGeometry::Point* points3[25];
  for (int i = 0; i < num3; ++i) {
    points3[i] = new Amanzi::AmanziGeometry::Point(x, y, z);
  }

  // Now delete some points in 3D - point 3, 4, 5, 8, 7, 16, 18, 17
  // Deletion of 3 is a stand alone deletion
  // Deletion of 4 followed by point 5 will cause merging of two holes
  // Deletion of 8 followed by point 7 will cause merging of two holes
  // Deletion of 16, 18, 17 will cause bridging of three holes

  // no obvious way of checking the underlying coordinate buffer
  // machinery except through a debugger

  delete points3[2];

  delete points3[4];
  delete points3[5];

  delete points3[8];
  delete points3[7];

  delete points3[16];
  delete points3[18];
  delete points3[17];


  // Create a bunch of mixed dim points

  int num = 25;
  std::vector<Amanzi::AmanziGeometry::Point*> points;
  points.resize(25);

  for (int i = 0; i < num; ++i) {
    if (i % 2)
      points[i] = new Amanzi::AmanziGeometry::Point(x, y);
    else {
      points[i] = new Amanzi::AmanziGeometry::Point(3);
      points[i]->set(x, y, z);
    }
    x = x + 1.0;
    y = y + 0.5;
    z = z + 2.0;
  }


  // Don't bother deleting other points. They will get deleted when exiting


  // we should test the mathematical operations as well here
}
