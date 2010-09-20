#ifndef _ENTITY_MAPS_H_
#define _ENTITY_MAPS_H_

#include "dbc.hh"
#include "Entity_kind.hh"

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <map>


namespace STK_mesh
{
/* Class Entity_map
 *
 * Models the relationship between STK_mesh entity labels and
 * Mesh_data labels and a consecutive indexed set of Mesh_data labels
 * (to facilitate looping over entity kinds).
 */

class Entity_map
{

    unsigned int dimension_;

    stk::mesh::EntityRank ranks_ [3];
    Mesh_data::Entity_kind kinds_ [3];

    std::map<stk::mesh::EntityRank, Mesh_data::Entity_kind> rank_to_kind_;
    std::map<Mesh_data::Entity_kind, stk::mesh::EntityRank> kind_to_rank_;

    static bool valid_dimension_ (unsigned int dimension);

    inline void create_maps_ ();

public:

    Entity_map (int dimension) : dimension_ (dimension)
    {
        ASSERT (valid_dimension_ (dimension));
        create_maps_ ();
    }

    unsigned int dimension () const { return dimension_; }

    stk::mesh::EntityRank  kind_to_rank (Mesh_data::Entity_kind kind) const;
    Mesh_data::Entity_kind rank_to_kind (stk::mesh::EntityRank   rank) const;

};

void Entity_map::create_maps_ ()
{

    rank_to_kind_ [stk::mesh::Node] = Mesh_data::NODE;
    rank_to_kind_ [stk::mesh::Edge] = (dimension_ == 2) ? Mesh_data::FACE : Mesh_data::EDGE;
    rank_to_kind_ [stk::mesh::Face] = (dimension_ == 2) ? Mesh_data::CELL : Mesh_data::FACE;
    rank_to_kind_ [stk::mesh::Element] = Mesh_data::CELL;


    kind_to_rank_ [Mesh_data::NODE] = stk::mesh::Node;
    kind_to_rank_ [Mesh_data::EDGE] = stk::mesh::Edge;
    kind_to_rank_ [Mesh_data::FACE] = (dimension_ == 2) ? stk::mesh::Edge : stk::mesh::Face;
    kind_to_rank_ [Mesh_data::CELL] = (dimension_ == 2) ? stk::mesh::Face : stk::mesh::Element;

}


}




#endif /* _ENTITY_MAPS_H_ */
