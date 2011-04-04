#ifndef _CELL_TOPOLOGY_HH_
#define _CELL_TOPOLOGY_HH_


#include "Element_types.hh"
#include "Element_field_types.hh"

namespace shards { struct CellTopology; }

namespace STK_mesh
{

const shards::CellTopology get_topology_data (const Mesh_data::ELEMENT_TYPE&);

}
#endif
