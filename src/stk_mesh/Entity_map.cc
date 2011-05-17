/* Copyright 2011 LANS */

#include "Entity_map.hh"

namespace Amanzi {
namespace AmanziMesh {
namespace STK {

bool Entity_map::valid_dimension_(unsigned int dimension) {
    return (dimension == 2 || dimension == 3);
}

stk::mesh::EntityRank
Entity_map::kind_to_rank(Entity_kind kind) const {
    ASSERT(entity_valid_kind(kind));
    return kind_to_rank_.find (kind)->second;
}


Entity_kind
Entity_map::rank_to_kind(stk::mesh::EntityRank rank) const {
    return rank_to_kind_.find (rank)->second;
}

} // close namespace STK 
} // close namespace Mesh 
} // close namespace Amanzi 

