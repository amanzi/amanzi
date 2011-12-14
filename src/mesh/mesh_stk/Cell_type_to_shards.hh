#ifndef _CELL_TOPOLOGY_HH_
#define _CELL_TOPOLOGY_HH_


#include "MeshDefs.hh"
#include "Element_field_types.hh"

namespace shards { 
struct CellTopology; 
}

namespace Amanzi {
namespace AmanziMesh {
namespace STK {

const shards::CellTopology get_topology_data (const Cell_type& t);

} // close namespace STK 
} // close namespace Mesh 
} // close namespace Amanzi 

#endif
