#include <iostream>
#include <boost/format.hpp>
#include <boost/lambda/lambda.hpp>
namespace bl = boost::lambda;

#include "Mesh_STK_Impl.hh"
#include "dbc.hh"
#include "Mesh_common.hh"
#include "GeometricModel.hh"
#include "LabeledSetRegion.hh"


#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/SetOwners.hpp>

namespace Amanzi {
namespace AmanziMesh {
namespace STK {

// Constructors
// ------------

Mesh_STK_Impl::Mesh_STK_Impl (int space_dimension, 
                              const Epetra_MpiComm *comm,
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

unsigned int Mesh_STK_Impl::num_sets (stk::mesh::EntityRank rank) const 
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
Mesh_STK_Impl::count_entities (stk::mesh::EntityRank rank, Parallel_type category) const
{
  Entity_vector e;
  get_entities(rank, category, e);
  return e.size();
}


unsigned int 
Mesh_STK_Impl::count_entities (const stk::mesh::Part& part, Parallel_type category) const
{
  Entity_vector e;
  get_entities(part, category, e);
  return e.size();
}

/** 
 * 
 * @param rank 
 * @param category 
 * @param entities 
 */
void 
Mesh_STK_Impl::get_entities (stk::mesh::EntityRank rank, Parallel_type category, Entity_vector& entities) const
{
  get_entities_ (selector_ (category), rank, entities);
}

void 
Mesh_STK_Impl::get_entities (const stk::mesh::Part& part, Parallel_type category, Entity_vector& entities) const
{
  const stk::mesh::Selector part_selector = part & selector_ (category);
  const stk::mesh::EntityRank rank  = part.primary_entity_rank ();
  get_entities_ (part_selector, rank, entities);
}

void 
Mesh_STK_Impl::get_entities_ (const stk::mesh::Selector& selector, stk::mesh::EntityRank rank,
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
Mesh_STK_Impl::element_to_faces (stk::mesh::EntityId element, Entity_Ids& ids) const
{
  // Look up element from global id.
  const int cell_rank = entity_map_->kind_to_rank (CELL);
  const int face_rank = entity_map_->kind_to_rank (FACE);

  stk::mesh::Entity *entity = id_to_entity(cell_rank, element);
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
Mesh_STK_Impl::element_to_face_dirs(stk::mesh::EntityId element, 
                                    std::vector<int>& dirs) const
{
  const int me = communicator_->MyPID();
  const int cell_rank = entity_map_->kind_to_rank (CELL);
  const int face_rank = entity_map_->kind_to_rank (FACE);

  stk::mesh::Entity *entity = id_to_entity(cell_rank, element);

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

void
Mesh_STK_Impl::element_to_faces_and_dirs(stk::mesh::EntityId element,
                                         Entity_Ids& ids,
                                         std::vector<int>& dirs) const
{
  // Look up element from global id.
  const int cell_rank = entity_map_->kind_to_rank (CELL);
  const int face_rank = entity_map_->kind_to_rank (FACE);

  stk::mesh::Entity *entity = id_to_entity(cell_rank, element);
  ASSERT (entity->identifier () == element);

  const CellTopologyData* topo = stk::mesh::get_cell_topology (*entity);

  ASSERT(topo != NULL);

  stk::mesh::PairIterRelation faces = entity->relations( face_rank );
  for (stk::mesh::PairIterRelation::iterator it = faces.begin (); it != faces.end (); ++it)
  {
    stk::mesh::EntityId gid(it->entity ()->identifier ());
    ids.push_back (gid);

    stk::mesh::FieldTraits<Id_field_type>::data_type *owner = 
        stk::mesh::field_data<Id_field_type>(*face_owner_, *(it->entity()));
    int dir(1);
    if (*owner != element) {
      dir = -1;
    }
    dirs.push_back (dir);    
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
Mesh_STK_Impl::element_to_nodes (stk::mesh::EntityId element, Entity_Ids& ids) const
{

  const int cell_rank = entity_map_->kind_to_rank (CELL);
  const int node_rank = entity_map_->kind_to_rank (NODE);

  stk::mesh::Entity *entity = id_to_entity(cell_rank, element);

  stk::mesh::PairIterRelation nodes = entity->relations (node_rank);

  for (stk::mesh::PairIterRelation::iterator it = nodes.begin (); it != nodes.end (); ++it)
  {
    ids.push_back (it->entity ()->identifier ());
  }
}

void Mesh_STK_Impl::face_to_nodes (stk::mesh::EntityId element, Entity_Ids& ids) const
{
  const int from_rank = entity_map_->kind_to_rank (FACE);
  const int to_rank = entity_map_->kind_to_rank (NODE);
  stk::mesh::Entity *entity = id_to_entity(from_rank, element);
    
  stk::mesh::PairIterRelation nodes = entity->relations (to_rank);
    
  for (stk::mesh::PairIterRelation::iterator it = nodes.begin (); it != nodes.end (); ++it)
  {
    ids.push_back (it->entity ()->identifier ());
  }
}

// -------------------------------------------------------------
// Mesh_STK_Impl::face_to_elements
// -------------------------------------------------------------
void
Mesh_STK_Impl::face_to_elements(stk::mesh::EntityId face, Entity_Ids& ids) const
{
  const int cell_rank = entity_map_->kind_to_rank (CELL);
  const int face_rank = entity_map_->kind_to_rank (FACE);

  stk::mesh::Entity *entity = id_to_entity(face_rank, face);
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
    
// -------------------------------------------------------------
// Mesh_STK_Impl::node_to_faces
// -------------------------------------------------------------
void
Mesh_STK_Impl::node_to_faces(stk::mesh::EntityId element, Entity_Ids& ids) const
{
  const int from_rank = entity_map_->kind_to_rank (NODE);
  const int to_rank = entity_map_->kind_to_rank (FACE);
  stk::mesh::Entity *entity = id_to_entity(from_rank, element);
    
  stk::mesh::PairIterRelation nodes = entity->relations (to_rank);
    
  for (stk::mesh::PairIterRelation::iterator it = nodes.begin (); it != nodes.end (); ++it)
  {
    ids.push_back (it->entity ()->identifier ());
  }
}

// -------------------------------------------------------------
// Mesh_STK_Impl::node_to_elements
// -------------------------------------------------------------
void
Mesh_STK_Impl::node_to_elements(stk::mesh::EntityId element, Entity_Ids& ids) const
{
  const int from_rank = entity_map_->kind_to_rank (NODE);
  const int to_rank = entity_map_->kind_to_rank (CELL);
  stk::mesh::Entity *entity = id_to_entity(from_rank, element);
    
  stk::mesh::PairIterRelation nodes = entity->relations (to_rank);
    
  for (stk::mesh::PairIterRelation::iterator it = nodes.begin (); it != nodes.end (); ++it)
  {
    ids.push_back (it->entity ()->identifier ());
  }
}



double const * Mesh_STK_Impl::coordinates (stk::mesh::Entity* node) const
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
Mesh_STK_Impl::coordinates (stk::mesh::EntityId node) const
{

  stk::mesh::Entity *entity = id_to_entity(stk::mesh::Node, node);
  return coordinates (entity);

}

void
Mesh_STK_Impl::set_coordinates (stk::mesh::EntityId node, const double *coords) 
{

  stk::mesh::Entity *entity = id_to_entity(stk::mesh::Node, node);
  double * node_coordinates = stk::mesh::field_data(coordinate_field_, *entity);

  std::copy(coords, coords+space_dimension_, node_coordinates);
}

// -------------------------------------------------------------
// Mesh_STK_Impl::id_to_entity
// -------------------------------------------------------------
stk::mesh::Entity* 
Mesh_STK_Impl::id_to_entity (stk::mesh::EntityRank rank, 
                             stk::mesh::EntityId id) const
{
  stk::mesh::Entity *entity = bulk_data_->get_entity(rank, id);

  ASSERT(entity != NULL);
  ASSERT(entity->identifier() == id);
  return entity;
}


stk::mesh::Part* 
Mesh_STK_Impl::get_set (unsigned int set_id, stk::mesh::EntityRank rank)
{
    
  Id_map::const_iterator part_it = set_to_part_.find (std::make_pair (rank, set_id));
  ASSERT (part_it != set_to_part_.end ());

  return part_it->second;
  
}

stk::mesh::Part* Mesh_STK_Impl::get_set (const char* name, stk::mesh::EntityRank rank)
{
  stk::mesh::Part *part = meta_data_->get_part (name);
  if (part)
    ASSERT (part->primary_entity_rank () == rank);

  return part;
}

stk::mesh::Part* Mesh_STK_Impl::get_set (const std::string name, stk::mesh::EntityRank rank)
{

  stk::mesh::Part *part = meta_data_->get_part (name);
  if (part)
    ASSERT (part->primary_entity_rank () == rank);

  return part;
}

void Mesh_STK_Impl::get_sets (stk::mesh::EntityRank rank, stk::mesh::PartVector& sets) const
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

void Mesh_STK_Impl::get_set_ids (stk::mesh::EntityRank rank, std::vector<unsigned int> &ids) const
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

stk::mesh::Selector Mesh_STK_Impl::selector_ (Parallel_type category) const
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

bool Mesh_STK_Impl::valid_id (unsigned int id, stk::mesh::EntityRank rank) const
{
  return (set_to_part_.find (std::make_pair (rank, id)) != set_to_part_.end ());
}



// Static Information & Validators
// -------------------------------

stk::mesh::EntityRank Mesh_STK_Impl::get_element_type (int space_dimension)
{
  ASSERT (valid_dimension (space_dimension));
  return (space_dimension == 3) ? stk::mesh::Element : stk::mesh::Face;
}

stk::mesh::EntityRank Mesh_STK_Impl::get_face_type (int space_dimension)
{
  ASSERT (valid_dimension (space_dimension));
  return (space_dimension == 3) ? stk::mesh::Face : stk::mesh::Edge;
}

bool Mesh_STK_Impl::valid_dimension (int space_dimension)
{
  return (space_dimension >= 2) && (space_dimension <= 3);
}

bool Mesh_STK_Impl::valid_rank (stk::mesh::EntityRank rank)
{
  return (rank == stk::mesh::Node || rank == stk::mesh::Edge ||
          rank == stk::mesh::Face || rank == stk::mesh::Element);
}


void Mesh_STK_Impl::add_set_part_relation_ (unsigned int set_id, stk::mesh::Part& part)
{
  const unsigned int part_id = part.mesh_meta_data_ordinal ();
  const stk::mesh::EntityRank rank = part.primary_entity_rank ();
  const Rank_and_id rank_set_id = std::make_pair (rank, set_id);

  ASSERT (set_to_part_.find (rank_set_id) == set_to_part_.end ());

  set_to_part_ [rank_set_id]  = &part;

}


// Object validators
// -----------------


bool Mesh_STK_Impl::dimension_ok_ () const
{
  return valid_dimension (space_dimension_);
}


// -------------------------------------------------------------
// Mesh_STK_Impl::summary
// -------------------------------------------------------------
void
Mesh_STK_Impl::summary(std::ostream& os) const
{
  const int nproc(communicator_->NumProc());
  const int me(communicator_->MyPID());

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
    communicator_->Barrier();
  }
  
}

// -------------------------------------------------------------
// Mesh_STK_Impl::redistribute
// -------------------------------------------------------------
/** 
 * This redistributes ownership of cells in the mesh according to the
 * specified @c cellmap.  
 * 
 * @param cellmap @em 1-based map of cell indexes
 */
void
Mesh_STK_Impl::redistribute(const Epetra_Map& cellmap)
{
  stk::mesh::Selector owned(meta_data_->locally_owned_part());
  stk::mesh::EntityVector cells;
  stk::mesh::get_selected_entities(owned, bulk_data_->buckets(stk::mesh::Element), cells);

  std::vector<int> gid(cells.size());

  int i(0);
  stk::mesh::EntityVector::iterator c;
  for (c = cells.begin(), i = 0; c != cells.end(); c++, i++) {
    gid[i] = (*c)->identifier();
  }

  std::vector<int> pid(cells.size());
  std::vector<int> lid(cells.size());

  cellmap.RemoteIDList(cells.size(), &gid[0], &pid[0], &lid[0]);

  std::vector<stk::mesh::EntityProc> eproc;

  for (c = cells.begin(), i = 0; c != cells.end(); c++, i++) {
    eproc.push_back(stk::mesh::EntityProc(*c, pid[i]));
  }

  // move the cells around

  bulk_data_->modification_begin();
  bulk_data_->change_entity_owner(eproc);
  bulk_data_->modification_end();

  // move the other stuff around
    
  bulk_data_->modification_begin();
  stk::mesh::set_owners<stk::mesh::LowestRankSharingProcOwns>(*bulk_data_);
  bulk_data_->modification_end();

  stk::mesh::check_face_ownership(*bulk_data_);
  stk::mesh::check_node_ownership(*bulk_data_);

    
}

} // close namespace STK 
} // close namespace Mesh 
} // close namespace Amanzi 

