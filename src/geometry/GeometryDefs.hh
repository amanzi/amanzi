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

typedef int Entity_ID; // should be consistent with similar definition in class Mesh
using Point_List = std::vector<AmanziGeometry::Point>;
template <typename T>
using View_type = std::vector<T>;
using Point_View = View_type<AmanziGeometry::Point>;
using Entity_ID_View = View_type<Entity_ID>;

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

std::string inline to_string(const RegionType rtype)
{
  switch (rtype) {
  case (RegionType::BOX): {
    return "region: box";
  } break;
  case (RegionType::PLANE): {
    return "region: plane";
  } break;
  case (RegionType::LABELEDSET): {
    return "region: labeled set";
  } break;
  case (RegionType::LAYER): {
    return "region: layer";
  } break;
  case (RegionType::SURFACE): {
    return "region: surface";
  } break;
  case (RegionType::POINT): {
    return "region: point";
  } break;
  case (RegionType::COLORFUNCTION): {
    return "region: color function";
  } break;
  case (RegionType::LOGICAL): {
    return "region: logical";
  } break;
  case (RegionType::POLYGON): {
    return "region: polygon";
  } break;
  case (RegionType::ENUMERATED): {
    return "region: enumerated";
  } break;
  case (RegionType::BOUNDARY): {
    return "region: boundary";
  } break;
  case (RegionType::BOX_VOF): {
    return "region: box vof";
  } break;
  case (RegionType::LINE_SEGMENT): {
    return "region: line segment";
  } break;
  case (RegionType::CYLINDER): {
    return "region: cylinder";
  } break;
  case (RegionType::ALL): {
    return "region: all";
  } break;
  default: {
    return "unknown";
  }
  }
}

inline std::ostream&
operator<<(std::ostream& os, const RegionType& rtype)
{
  os << to_string(rtype);
  return os;
}


enum class LifeCycleType { PERMANENT = 0, TEMPORARY };


enum class BoolOpType { NOBOOLEAN = -1, COMPLEMENT, UNION, INTERSECT, SUBTRACT };

const double TOL = 1.0e-08;

// arbitrary number to avoid clashing
// with IDs of LabeledSet regions
const unsigned int REGION_ID_OFFSET = 59049;

} // namespace AmanziGeometry
} // namespace Amanzi

#endif
