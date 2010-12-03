#include "Mesh.hh"
#include "dbc.hh"

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>

namespace STK_mesh
{

// Constructors
// ------------

Mesh::Mesh (int space_dimension, 
            const Epetra_MpiComm& comm,
            Entity_map* entity_map,
            stk::mesh::MetaData *meta_data,
            stk::mesh::BulkData *bulk_data,
            const Id_map& set_to_part,
            Vector_field_type &coordinate_field) :
    space_dimension_ (space_dimension),
    communicator_ (comm),
    entity_map_ (entity_map),
    meta_data_ (meta_data),
    bulk_data_ (bulk_data),
    set_to_part_ (set_to_part),
    coordinate_field_ (coordinate_field)
{

    ASSERT (dimension_ok_ ());
    ASSERT (meta_data_.get ());
    ASSERT (bulk_data_.get ());
}

// Information Getters
// -------------------

unsigned int Mesh::num_sets (stk::mesh::EntityRank rank) const 
{
    int count = 0;

    for (Id_map::const_iterator it = set_to_part_.begin ();
         it != set_to_part_.end ();
         ++it)
    {
        if (it->second->primary_entity_rank () == rank) ++count;
    }

    return count;
}

/** 
 * Ideally, we'd like to trust stk::mesh to understand @c category ,
 * but it doesn't.  See get_entities() for how things really need to
 * be done.
 * 
 * @param rank 
 * @param category 
 * 
 * @return number of entities found
 */
unsigned int 
Mesh::count_entities (stk::mesh::EntityRank rank, Element_Category category) const
{
    Entity_vector e;
    get_entities(rank, category, e);
    return e.size();
}


unsigned int 
Mesh::count_entities (const stk::mesh::Part& part, Element_Category category) const
{
    Entity_vector e;
    get_entities(part, category, e);
    return e.size();
}

/** 
 * Because cell-to-face relations are problematic, we can't really
 * trust stk::mesh to get ghost cells and faces right.  
 * 
 * @param rank 
 * @param category 
 * @param entities 
 */
void 
Mesh::get_entities (stk::mesh::EntityRank rank, Element_Category category, Entity_vector& entities) const
{
    get_entities_ (selector_ (category), rank, entities);
}

void 
Mesh::get_entities (const stk::mesh::Part& part, Element_Category category, Entity_vector& entities) const
{
    const stk::mesh::Selector part_selector = part & selector_ (category);
    const stk::mesh::EntityRank rank  = part.primary_entity_rank ();
    get_entities_ (part_selector, rank, entities);
}

void 
Mesh::get_entities_ (const stk::mesh::Selector& selector, stk::mesh::EntityRank rank,
                     Entity_vector& entities) const
{
    stk::mesh::get_selected_entities (selector, bulk_data_->buckets (rank), entities);
}


/** 
 * This may only be safe if @c element is locally owned or shared.
 * 
 * @param element @em global element identifier (in)
 * @param ids @em global face identifiers (out)
 */
void 
Mesh::element_to_faces (stk::mesh::EntityId element, Entity_Ids& ids) const
{
    // Look up element from global id.
    const int node_rank = entity_map_->kind_to_rank (Mesh_data::NODE);
    const int cell_rank = entity_map_->kind_to_rank (Mesh_data::CELL);
    const int face_rank = entity_map_->kind_to_rank (Mesh_data::FACE);

    stk::mesh::Entity *entity = bulk_data_->get_entity(cell_rank, element);
    ASSERT (entity->identifier () == element);

    const CellTopologyData* topo = stk::mesh::get_cell_topology (*entity);

    ASSERT(topo != NULL);

    stk::mesh::PairIterRelation faces = entity->relations( face_rank );
    for (stk::mesh::PairIterRelation::iterator it = faces.begin (); it != faces.end (); ++it)
    {
        ids.push_back (it->entity ()->identifier ());
    }

    ASSERT (ids.size () == topo->side_count);
}

/** 
 * This may only be safe if @c element is locally owned or shared.
 * 
 * @param element @em global element identifier (in)
 * @param ids @em global node identifiers (out)
 */
void 
Mesh::element_to_nodes (stk::mesh::EntityId element, Entity_Ids& ids) const
{

    const int cell_rank = entity_map_->kind_to_rank (Mesh_data::CELL);
    const int node_rank = entity_map_->kind_to_rank (Mesh_data::NODE);

    stk::mesh::Entity *entity = bulk_data_->get_entity(cell_rank, element);

    stk::mesh::PairIterRelation nodes = entity->relations (node_rank);

    for (stk::mesh::PairIterRelation::iterator it = nodes.begin (); it != nodes.end (); ++it)
    {
        ids.push_back (it->entity ()->identifier ());
    }

    ASSERT (ids.size () == 8);

}

void Mesh::face_to_nodes (stk::mesh::EntityId element, Entity_Ids& ids) const
{
    const int from_rank = entity_map_->kind_to_rank (Mesh_data::FACE);
    const int to_rank = entity_map_->kind_to_rank (Mesh_data::NODE);
    stk::mesh::Entity *entity = bulk_data_->get_entity(from_rank, element);
    
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

/** 
 * 
 * 
 * @param node @em global node identifier
 * 
 * @return node coordinates
 */
double const * 
Mesh::coordinates (stk::mesh::EntityId node) const
{

    stk::mesh::Entity *entity = bulk_data_->get_entity(stk::mesh::Node, node);
    return coordinates (entity);

}

stk::mesh::Part* 
Mesh::get_set (unsigned int set_id, stk::mesh::EntityRank rank)
{
    
    Id_map::const_iterator part_it = set_to_part_.find (std::make_pair (rank, set_id));
    ASSERT (part_it != set_to_part_.end ());

    return part_it->second;

}

stk::mesh::Part* Mesh::get_set (const char* name, stk::mesh::EntityRank rank)
{

    stk::mesh::Part *part = meta_data_->get_part (name);
    ASSERT (part);
    ASSERT (part->primary_entity_rank () == rank);

    return part;
}

void Mesh::get_sets (stk::mesh::EntityRank rank, stk::mesh::PartVector& sets) const
{
    
    ASSERT (sets.size () == 0);

    for (Id_map::const_iterator it = set_to_part_.begin ();
         it != set_to_part_.end ();
         ++it)
    {
        if (it->first.first == rank) sets.push_back (it->second);
    }

    ASSERT (sets.size () == num_sets (rank));

}

void Mesh::get_set_ids (stk::mesh::EntityRank rank, std::vector<unsigned int> &ids) const
{
    ASSERT (ids.size () == 0);
    
    for (Id_map::const_iterator it = set_to_part_.begin ();
         it != set_to_part_.end ();
         ++it)
    {
        if (it->first.first == rank) ids.push_back (it->first.second);
    }

}


// Manipulators
// ------------

stk::mesh::Selector Mesh::selector_ (Element_Category category) const
{
    ASSERT (valid_category (category));

    stk::mesh::Selector s;
    switch (category) {
    case (OWNED):
      s |= meta_data_->locally_owned_part();
      break;
    case (GHOST):
      s |= meta_data_->globally_shared_part();
      break;
    case (USED):
      s |= meta_data_->locally_owned_part();
      s |= meta_data_->globally_shared_part();
      break;
    }
    return s;
}



// Argument validators
// -------------------

bool Mesh::valid_id (unsigned int id, stk::mesh::EntityRank rank) const
{
    return (set_to_part_.find (std::make_pair (rank, id)) != set_to_part_.end ());
}



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

