#include "Cell_topology.hh"

#include "Element_types.hh"
#include "Element_field_types.hh"

#include "Shards_CellTopology.hpp"

#include "dbc.hh"

#include "stk_mesh_error.hh"

namespace STK_mesh
{


const shards::CellTopology get_topology_data (const Mesh_data::ELEMENT_TYPE& type)
{
    ASSERT (Mesh_data::ok_type (type));

    switch (type) {
    case (Mesh_data::HEX):
        return shards::getCellTopologyData< shards::Hexahedron<8> >();
        break;
    case (Mesh_data::TETRA):
        return shards::getCellTopologyData< shards::Tetrahedron<4> >();
        break;
    case (Mesh_data::PYRAMID):
        return shards::getCellTopologyData< shards::Pyramid<5> >();
        break;
    case (Mesh_data::WEDGE):
        return shards::getCellTopologyData< shards::Wedge<6> >();
    default:
        // fall through
        break;
    }

    Exceptions::amanzi_throw( STKMeshError ("Element not supported yet") );

}


}
