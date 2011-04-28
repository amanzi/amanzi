/* Copyright 2011 LANS */

#include "Entity_map.hh"

namespace STK_mesh {

bool Entity_map::valid_dimension_(unsigned int dimension) {
    return (dimension == 2 || dimension == 3);
}

stk::mesh::EntityRank
Entity_map::kind_to_rank(Mesh_data::Entity_kind kind) const {
    ASSERT(Mesh_data::valid_entity_kind(kind));
    return kind_to_rank_.find (kind)->second;
}


Mesh_data::Entity_kind
Entity_map::rank_to_kind(stk::mesh::EntityRank rank) const {
    return rank_to_kind_.find (rank)->second;
}
}
