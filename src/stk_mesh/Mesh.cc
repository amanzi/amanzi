#include <iostream>
#include <boost/format.hpp>
#include <boost/lambda/lambda.hpp>
namespace bl = boost::lambda;

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
    coordinate_field_ (coordinate_field),
    face_owner_(meta_data_->get_field<Id_field_type>("FaceOwner"))
{

    ASSERT (dimension_ok_ ());
    ASSERT (meta_data_.get ());
    ASSERT (bulk_data_.get ());
    ASSERT (face_owner_ != NULL);
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

// -------------------------------------------------------------
// Mesh::remove_bad_ghosts_
// -------------------------------------------------------------
/** 
 * This routine takes a list of entities and removes those ghost
 * entities that are not appropriate.  
 *
 * Nodes: stk::mesh ghosts all nodes in a ghosted cell, even if you
 * tell it not to.  So, a node must be included in a cell owned by
 * this process.  Others are removed.
 *
 * Faces: stk::mesh ghosts all faces in a ghosted cell, even if you
 * tell it not to.  So each remotely owned face is checked to see if
 * it's related to a locally-owned cell.  If not, it's removed from
 * the list.
 *
 * FIXME: this appears to be unnecessary; it was based on my
 * misconception of what needed to be ghosted.
 * 
 * @param entities 
 */
void
Mesh::remove_bad_ghosts_(Entity_vector& entities) const
{
    const int me(communicator_.MyPID());
    Entity_vector tmp(entities);

    for (Entity_vector::iterator i = tmp.begin(); i != tmp.end(); i++) {
        bool bad(false);
        if ((*i)->owner_rank() != me) {

            if ((*i)->entity_rank() == stk::mesh::Node) {
                stk::mesh::EntityVector thenode, cells;
                thenode.push_back(*i);
                stk::mesh::get_entities_through_relations(thenode, stk::mesh::Element, cells);
                
                for (stk::mesh::EntityVector:: iterator c = cells.begin(); c != cells.end(); c++) {
                    if ((*c)->owner_rank() == me) {
                        bad = false;
                        break;
                    } else {
                        bad = true;
                    }
                }
            }

            if ((*i)->entity_rank() == stk::mesh::Face) { 
                stk::mesh::EntityVector theface, cells;
                theface.push_back(*i);
                stk::mesh::get_entities_through_relations(theface, stk::mesh::Element, cells);
                
                
                // ghost faces must relate to more than one cell
                bad = bad || (cells.size() < 2);
                
                // ghost faces must relate to a cell owned by this process
                bad = bad ||  (cells.front()->owner_rank() != me &&
                               cells.back()->owner_rank() != me); 
                
                // std::cerr << "Process " << me << ": Face " << (*i)->identifier() 
                //           << " connects cell " << cells.front()->identifier() 
                //           << " to "  << cells.back()->identifier() << " " << bad << std::endl;
                
            }
        }
        if (bad) entities.erase(std::remove(entities.begin(), entities.end(), *i));
    }
}


/** 
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
    // remove_bad_ghosts_(entities);
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
    const int cell_rank = entity_map_->kind_to_rank (Mesh_data::CELL);
    const int face_rank = entity_map_->kind_to_rank (Mesh_data::FACE);

    stk::mesh::Entity *entity = bulk_data_->get_entity(cell_rank, element);
    ASSERT (entity->identifier () == element);

    const CellTopologyData* topo = stk::mesh::get_cell_topology (*entity);

    ASSERT(topo != NULL);

    stk::mesh::PairIterRelation faces = entity->relations( face_rank );
    for (stk::mesh::PairIterRelation::iterator it = faces.begin (); it != faces.end (); ++it)
    {
        stk::mesh::EntityId gid(it->entity ()->identifier ());
        ids.push_back (gid);
    }

    ASSERT (ids.size () == topo->side_count);
}

void
Mesh::element_to_face_dirs(stk::mesh::EntityId element, 
                           std::vector<int>& dirs) const
{
    const int me = communicator_.MyPID();
    const int cell_rank = entity_map_->kind_to_rank (Mesh_data::CELL);
    const int face_rank = entity_map_->kind_to_rank (Mesh_data::FACE);

    stk::mesh::Entity *entity = bulk_data_->get_entity(cell_rank, element);

    ASSERT (entity->identifier () == element);

    stk::mesh::PairIterRelation faces = entity->relations( face_rank );
    for (stk::mesh::PairIterRelation::iterator it = faces.begin (); it != faces.end (); ++it)
    {
        stk::mesh::FieldTraits<Id_field_type>::data_type *owner = 
            stk::mesh::field_data<Id_field_type>(*face_owner_, *(it->entity()));
        int dir(1);
        if (*owner != element) {
            dir = -1;
        }
        dirs.push_back (dir);
    }
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

// -------------------------------------------------------------
// Mesh::face_to_elements
// -------------------------------------------------------------
void
Mesh::face_to_elements(stk::mesh::EntityId face, Entity_Ids& ids) const
{
    const int cell_rank = entity_map_->kind_to_rank (Mesh_data::CELL);
    const int face_rank = entity_map_->kind_to_rank (Mesh_data::FACE);

    stk::mesh::Entity *entity = bulk_data_->get_entity(face_rank, face);
    ASSERT (entity->identifier () == face);

    stk::mesh::PairIterRelation cells = entity->relations( cell_rank );
    ASSERT(!cells.empty());

    for (stk::mesh::PairIterRelation::iterator it = cells.begin (); it != cells.end (); ++it)
    {
        stk::mesh::EntityId gid(it->entity ()->identifier ());
        ids.push_back (gid);
    }

    ASSERT(ids.size() <= 2);
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
    ASSERT (category >= OWNED && category <= USED);

    stk::mesh::Selector owned(meta_data_->locally_owned_part());
    stk::mesh::Selector s;
    switch (category) {
    case (OWNED):
      s = owned;
      break;
    case (GHOST):
      s = !owned;
      // s &= meta_data_->globally_shared_part();
      break;
    case (USED):
      s = meta_data_->universal_part();
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


// -------------------------------------------------------------
// Mesh::summary
// -------------------------------------------------------------
void
Mesh::summary(std::ostream& os) const
{
    const int nproc(communicator_.NumProc());
    const int me(communicator_.MyPID());

    for (int p = 0; p < nproc; p++) {
        if (p == me) {
            os << boost::str(boost::format("Process %d: Nodes: %5d owned, %5d ghost, %5d used") %
                             me % 
                             count_entities(stk::mesh::Node, OWNED) %
                             count_entities(stk::mesh::Node, GHOST) %
                             count_entities(stk::mesh::Node, USED))
               << std::endl;

            os << boost::str(boost::format("Process %d: Faces: %5d owned, %5d ghost, %5d used") %
                             me % 
                             count_entities(stk::mesh::Face, OWNED) %
                             count_entities(stk::mesh::Face, GHOST) %
                             count_entities(stk::mesh::Face, USED))
               << std::endl;
            os << boost::str(boost::format("Process %d: Cells: %5d owned, %5d ghost, %5d used") %
                             me % 
                             count_entities(stk::mesh::Element, OWNED) %
                             count_entities(stk::mesh::Element, GHOST) %
                             count_entities(stk::mesh::Element, USED))
               << std::endl;
        }
        communicator_.Barrier();
    }
  
}

} // close namespace STK_mesh

