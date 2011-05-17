#ifndef _ENTITY_MAPS_H_
#define _ENTITY_MAPS_H_

#include "dbc.hh"
#include "MeshDefs.hh"

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <map>


namespace Amanzi {
namespace AmanziMesh {
namespace STK {

/* Class Entity_map
 *
 * Models the relationship between STK_mesh entity labels and
 * Mesh_data labels and a consecutive indexed set of Mesh_data labels
 * (to facilitate looping over entity kinds).
 */

class Entity_map
{

    unsigned int dimension_;

    stk::mesh::EntityRank ranks_  [3];
    Entity_kind kinds_ [3];

    std::map<stk::mesh::EntityRank, Entity_kind> rank_to_kind_;
    std::map<Entity_kind, stk::mesh::EntityRank> kind_to_rank_;

    static bool valid_dimension_ (unsigned int dimension);

    inline void create_maps_ ();

public:

    Entity_map (int dimension) : dimension_ (dimension)
    {
        ASSERT (valid_dimension_ (dimension));
        create_maps_ ();
    }

    unsigned int dimension () const { return dimension_; }

    stk::mesh::EntityRank  kind_to_rank (Entity_kind kind) const;
    Entity_kind rank_to_kind (stk::mesh::EntityRank   rank) const;

};

void Entity_map::create_maps_ ()
{

    rank_to_kind_ [stk::mesh::Node] = NODE;
    rank_to_kind_ [stk::mesh::Edge] = (dimension_ == 2) ? FACE : EDGE;
    rank_to_kind_ [stk::mesh::Face] = (dimension_ == 2) ? CELL : FACE;
    rank_to_kind_ [stk::mesh::Element] = CELL;


    kind_to_rank_ [NODE] = stk::mesh::Node;
    kind_to_rank_ [EDGE] = stk::mesh::Edge;
    kind_to_rank_ [FACE] = (dimension_ == 2) ? stk::mesh::Edge : stk::mesh::Face;
    kind_to_rank_ [CELL] = (dimension_ == 2) ? stk::mesh::Face : stk::mesh::Element;

}


} // close namespace STK 
} // close namespace Mesh 
} // close namespace Amanzi 




#endif /* _ENTITY_MAPS_H_ */
