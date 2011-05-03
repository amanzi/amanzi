/* Copyright 2011 LANS */

#include "Cell_topology.hh"

#include "Element_types.hh"
#include "Element_field_types.hh"

#include "Shards_CellTopology.hpp"

#include "dbc.hh"

namespace STK_mesh {

const shards::CellTopology
get_topology_data(const Mesh_data::ELEMENT_TYPE& type) {
    ASSERT(Mesh_data::ok_type(type));

    if (Mesh_data::TRIANGLE == type) {
        return shards::getCellTopologyData<shards::Triangle<3> >();
    } else if (Mesh_data::HEX == type) {
        return shards::getCellTopologyData<shards::Hexahedron<8> >();
    } else if (Mesh_data::TETRA == type) {
        return shards::getCellTopologyData<shards::Tetrahedron<4> >();
    }

    throw "Element not supported yet";
}
}
