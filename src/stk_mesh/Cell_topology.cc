/* Copyright 2011 LANS */

#include "Cell_topology.hh"

#include "Element_types.hh"
#include "Element_field_types.hh"
#include "stk_mesh_error.hh"
#include "Shards_CellTopology.hpp"

#include "dbc.hh"

namespace Amanzi {
namespace AmanziMesh {
namespace STK {


const shards::CellTopology
get_topology_data(const Cell_type& type) {
  ASSERT (cell_valid_type (type));

  switch (type) {
    case (HEX):
      return shards::getCellTopologyData< shards::Hexahedron<8> >();
      break;
    case (TET):
      return shards::getCellTopologyData< shards::Tetrahedron<4> >();
      break;
    case (PYRAMID):
      return shards::getCellTopologyData< shards::Pyramid<5> >();
      break;
    case (PRISM):
      return shards::getCellTopologyData< shards::Wedge<6> >();
      break;
    default:
      // fall through
      break;
  }
  Exceptions::amanzi_throw( STK::Error ("Element not supported yet") );
}

} // close namespace STK 
} // close namespace Mesh 
} // close namespace Amanzi 

