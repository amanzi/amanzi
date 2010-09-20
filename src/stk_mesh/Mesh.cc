#include "Mesh.hh"
#include "dbc.hh"

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>

// Helper functions and classes
namespace {

struct Element_has_id
{
    stk::mesh::EntityId id_;
    Element_has_id (stk::mesh::EntityId id) : id_ (id) { }
    bool operator () (stk::mesh::Entity *entity) { return (entity->identifier () == id_); }
};

}

namespace STK_mesh
{

// Constructors
// ------------

Mesh::Mesh (int space_dimension, 
            const Epetra_MpiComm& comm,
            Entity_map* entity_map,
            stk::mesh::MetaData *meta_data,
            stk::mesh::BulkData *bulk_data,
            Vector_field_type &coordinate_field) :
    space_dimension_ (space_dimension),
    communicator_ (comm),
    entity_map_ (entity_map),
    meta_data_ (meta_data),
    bulk_data_ (bulk_data),
    coordinate_field_ (coordinate_field)
{

    ASSERT (dimension_ok_ ());
    ASSERT (meta_data_.get ());
    ASSERT (bulk_data_.get ());

    update_ ();
}

// Information Getters
// -------------------

const stk::mesh::Selector& Mesh::selector_ (Element_Category category) const
{
    ASSERT (valid_category (category));

    if (category == OWNED) return owned_selector_;
    if (category == USED)  return used_selector_;
    if (category == GHOST) return ghost_selector_;

    throw "Invalid element category in Mesh::selector_";
}

stk::mesh::Entity* Mesh::id_to_entity (stk::mesh::EntityRank rank, 
                                       stk::mesh::EntityId id, 
                                       Element_Category category) const
{
    ASSERT (valid_rank (rank));
    ASSERT (valid_category (category));

    // Get the elements of correct rank
    Entity_vector entities;
    get_entities (rank, category, entities);

    // Search for the given id.
    Entity_vector::const_iterator entity = 
        std::find_if (entities.begin (), entities.end (), Element_has_id (id));

    if (entity == entities.end ()) throw "Unable to find specified entity in id_to_entity_";

    return *entity;
}

unsigned int Mesh::count_entities (stk::mesh::EntityRank rank, Element_Category category) const
{
    return stk::mesh::count_selected_entities (selector_ (category), bulk_data_->buckets (rank));
}

void Mesh::get_entities (stk::mesh::EntityRank rank, Element_Category category, Entity_vector& entities) const
{

    ASSERT (entities.size () == 0);
    stk::mesh::get_selected_entities (selector_ (category), bulk_data_->buckets (rank), entities);
    ASSERT (entities.size () == count_entities (rank, category));

}


void Mesh::element_to_faces (stk::mesh::EntityId element, Entity_Ids& ids) const
{
    // Look up element from global id.
    const int cell_rank = entity_map_->kind_to_rank (Mesh_data::CELL);
    const int face_rank = entity_map_->kind_to_rank (Mesh_data::FACE);
    stk::mesh::Entity *entity = id_to_entity (cell_rank, element, USED);
    ASSERT (entity->identifier () == element);

    // Get relation connections
    stk::mesh::PairIterRelation faces = entity->relations (face_rank);

    for (stk::mesh::PairIterRelation::iterator it = faces.begin (); it != faces.end (); ++it)
    {
        ids.push_back (it->entity ()->identifier ());
    }

    ASSERT (ids.size () == 6);
}

void Mesh::element_to_nodes (stk::mesh::EntityId element, Entity_Ids& ids) const
{

    const int cell_rank = entity_map_->kind_to_rank (Mesh_data::CELL);
    const int node_rank = entity_map_->kind_to_rank (Mesh_data::NODE);
    stk::mesh::Entity *entity = id_to_entity (cell_rank, element, USED);

    stk::mesh::PairIterRelation faces = entity->relations (node_rank);

    for (stk::mesh::PairIterRelation::iterator it = faces.begin (); it != faces.end (); ++it)
    {
        ids.push_back (it->entity ()->identifier ());
    }

    ASSERT (ids.size () == 8);

}

void Mesh::face_to_nodes (stk::mesh::EntityId element, Entity_Ids& ids) const
{
    const int from_rank = entity_map_->kind_to_rank (Mesh_data::FACE);
    const int to_rank = entity_map_->kind_to_rank (Mesh_data::NODE);
    stk::mesh::Entity *entity = id_to_entity (from_rank, element, USED);
    
    stk::mesh::PairIterRelation nodes = entity->relations (to_rank);
    
    for (stk::mesh::PairIterRelation::iterator it = nodes.begin (); it != nodes.end (); ++it)
    {
        ids.push_back (it->entity ()->identifier ());
    }

    ASSERT (ids.size () == 4);
    
}

double const * Mesh::coordinates (stk::mesh::Entity* node) const
{
    
    // Get an array of entity data.
    return field_data (coordinate_field_, *node);
}


double const * Mesh::coordinates (stk::mesh::EntityId node) const
{

    ASSERT (node > 0);
    ASSERT (node < count_entities (stk::mesh::Node, USED));

    stk::mesh::Entity *entity = id_to_entity (stk::mesh::Node, node, USED);
    return coordinates (entity);

}


// Manipulators
// ------------

void Mesh::update_ ()
{

    // Make sure whatever mesh data I'm caching is brought up-to-date here.

    // Update cached element selectors
    universal_selector_ = stk::mesh::Selector (meta_data_->universal_part ());
    owned_selector_     = stk::mesh::Selector (meta_data_->locally_owned_part ());
    ghost_selector_     = stk::mesh::Selector (meta_data_->globally_shared_part ());
    used_selector_      = owned_selector_ | ghost_selector_;
    
    notify_views_ ();

}


// void Mesh::rebalance_mesh (const Mesh::Entity_map& entity_map)
// {

//     // Loop over entities in entity_map
//     for (Entity_map::const_iterator map_it;
//          map_it != entity_map.end ();
//          ++map_it)
//     {

//         // Look up the element from the global id in the map
//         const stk::mesh::EntityId entity_global_id = map_it->first;
//         const unsigned int destination = map_it->second;

//         const stk::mesh::EntityId local_id = global_to_local (entity_global_id);
//         if (local_id != 0)
//         {
//             // If the element was found, check it's destination. If not on
//             // this rank, add it to the list of outgoing entities.
//             if (destination != rank_id ())
//             {

//             }

//         }        
        

//     }

// }

// Static Information & Validators
// -------------------------------

stk::mesh::EntityRank Mesh::get_element_type (int space_dimension)
{
    ASSERT (valid_dimension (space_dimension));
    return (space_dimension == 3) ? stk::mesh::Element : stk::mesh::Face;
}

stk::mesh::EntityRank Mesh::get_face_type (int space_dimension)
{
    ASSERT (valid_dimension (space_dimension));
    return (space_dimension == 3) ? stk::mesh::Face : stk::mesh::Edge;
}

bool Mesh::valid_dimension (int space_dimension)
{
    return (space_dimension >= 2) && (space_dimension <= 3);
}

bool Mesh::valid_rank (stk::mesh::EntityRank rank)
{
    return (rank == stk::mesh::Node || rank == stk::mesh::Edge ||
            rank == stk::mesh::Face || rank == stk::mesh::Element);
}



// Object validators
// -----------------


bool Mesh::dimension_ok_ () const
{
    return valid_dimension (space_dimension_);
}

}

