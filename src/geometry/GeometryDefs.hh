/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

// Definitions for geometry, regions.
//
// Author: Ethan Coon
//

#ifndef AMANZI_GEOMETRY_DEFS_HH_
#define AMANZI_GEOMETRY_DEFS_HH_

namespace Amanzi {
namespace AmanziGeometry {  

typedef int Entity_ID; 
typedef std::vector<Entity_ID> Entity_ID_List;

typedef int Set_ID;
typedef std::vector<Set_ID> Set_ID_List;

typedef enum {
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
  ALL
} RegionType;


typedef enum {
  PERMANENT=0,
  TEMPORARY
} LifeCycleType;


typedef enum {
  NOBOOLEAN=-1,
  COMPLEMENT,
  UNION,
  INTERSECT,
  SUBTRACT
} BoolOpType;

const double TOL = 1.0e-08;

// arbitrary number to avoid clashing
// with IDs of LabeledSet regions
const unsigned int REGION_ID_OFFSET = 59049;  


} // namespace
} // namespace




#endif
