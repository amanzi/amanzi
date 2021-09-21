/*
  Geometry

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Definitions for geometry, regions.
*/

#ifndef AMANZI_GEOMETRY_DEFS_HH_
#define AMANZI_GEOMETRY_DEFS_HH_

#include "Point.hh"

namespace Amanzi {
namespace AmanziGeometry {

typedef int Entity_ID;  // should be consistent with similar definition in class Mesh
typedef std::vector<Point> Point_List;

enum class RegionType {
  BOX,
  PLANE,
  LABELEDSET,
  LAYER,
  SURFACE,
  POINT,
  COLORFUNCTION,
  LOGICAL,
  POLYGON,
  ENUMERATED,
  BOUNDARY,
  BOX_VOF,
  LINE_SEGMENT,
  CYLINDER,
  ALL
};


enum class LifeCycleType {
  PERMANENT=0,
  TEMPORARY
};


enum class BoolOpType {
  NOBOOLEAN=-1,
  COMPLEMENT,
  UNION,
  INTERSECT,
  SUBTRACT
};

const double TOL = 1.0e-08;

// arbitrary number to avoid clashing
// with IDs of LabeledSet regions
const unsigned int REGION_ID_OFFSET = 59049;

}  // namespace AmanziGeometry
}  // namespace Amanzi

#endif
