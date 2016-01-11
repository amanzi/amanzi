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
typedef std::string Set_Name;



typedef enum {
  BOX=0,
  PLANE,
  LABELEDSET,
  LAYER,
  SURFACE,
  POINT,
  COLORFUNCTION,  
  LOGICAL,
  POLYGON,
  ENUMERATEDSET
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

} // namespace
} // namespace

#endif
