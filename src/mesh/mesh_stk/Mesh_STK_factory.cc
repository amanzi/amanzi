#include "dbc.hh"

#include <iostream>
#include <algorithm>
#include <iterator>
#include <sstream>

#include <boost/format.hpp>
#include <boost/lambda/lambda.hpp>
namespace bl = boost::lambda;

// Mesh_data
#include "Data.hh"
#include "Side_set.hh"
#include "Node_set.hh"
#include "Element_block.hh"
#include "Element_types.hh"
#include "Entity_map.hh"

// AmanziGeometry
#include "Geometry.hh"
#include "BoxRegion.hh"
#include "PointRegion.hh"

// STK_mesh
#include "Mesh_STK_factory.hh"
#include "Element_field_types.hh"
#include "Cell_type_to_shards.hh"
#include "Mesh_common.hh"
#include "stk_mesh_error.hh"

// Trilinos
#include <Shards_CellTopology.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologicalMetaData.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>

namespace Amanzi {
namespace AmanziMesh {
namespace STK {

Mesh_STK_factory::Mesh_STK_factory (const stk::ParallelMachine& comm, int bucket_size) 
    : parallel_machine_ (comm),
      bucket_size_ (bucket_size),
      communicator_ (comm)
{  }

/** 
 * This should only be called for serial 
 * 
 * @param data the mesh
 * @param fields data fields on the mesh
 * 
 * @return 
 */
Mesh_STK_Impl* Mesh_STK_factory::build_mesh (const Data::Data& data, 
					     const Data::Fields& fields,
					     const AmanziGeometry::GeometricModelPtr& gm)
{
  ASSERT(communicator_.NumProc() == 1);

  int ncell(data.parameters().num_elements_);
  Epetra_Map cmap(ncell, 1, communicator_);

  int nvert(data.parameters().num_nodes_);
  Epetra_Map vmap(nvert, 1, communicator_);

  return build_mesh(data, cmap, vmap, fields, gm);
}

/** 
 * Construct a Mesh instance, in parallel
 * 
 * @param data mesh description
 * @param cellmap map of local to global cell indexes
 * @param vertmap map of local to global vertex indexes
 * @param fields any data fields to put on the mesh
 * 
 * @return Mesh instance
 */
Mesh_STK_Impl* Mesh_STK_factory::build_mesh (const Data::Data& data, 
					     const Epetra_Map& cellmap,
					     const Epetra_Map& vertmap,
					     const Data::Fields& fields,
					     const AmanziGeometry::GeometricModelPtr& gm)
{

  // Update construction variables for the new mesh.
  const int space_dimension = data.parameters ().dimensions ();

  ASSERT (Mesh_STK_Impl::valid_dimension (space_dimension));

  entity_map_ = new Entity_map (space_dimension);
  meta_data_  = new stk::mesh::MetaData (stk::mesh::fem_entity_rank_names ());
  bulk_data_  = new stk::mesh::BulkData (*meta_data_, parallel_machine_, bucket_size_);

  node_rank_    = entity_map_->kind_to_rank (NODE);
  face_rank_    = entity_map_->kind_to_rank (FACE);
  element_rank_ = entity_map_->kind_to_rank (CELL);


  // Reset all of the mesh-specific data.
  face_id_ = 0;
  Parts (0).swap (element_blocks_);
  Parts (0).swap (side_sets_);
  Parts (0).swap (node_sets_);
  coordinate_field_ = 0;
  set_to_part_.clear ();
    
  // Build the data for the mesh object.

  build_meta_data_ (data, fields, gm);
  build_bulk_data_ (data, cellmap, vertmap, fields, gm);

  Mesh_STK_Impl *mesh(NULL);

  // NOTE: Ownership of the pointers meta_data_ and bulk_data_ are
  // transferred to the Mesh instance.  That's where the memory is
  // deleted. IMHO, smart pointers should be used ANYTIME memory is
  // allocated in one code block and deleted somewhere else. I
  // wasted too much of my life looking for memory leaks. -WAP

  mesh = new Mesh_STK_Impl (space_dimension, communicator_, entity_map_, 
                            meta_data_, bulk_data_,
                            set_to_part_,
                            *(meta_data_->get_field<Vector_field_type> (std::string ("coordinates"))));
    
  return mesh;
}


  void Mesh_STK_factory::build_meta_data_ (const Data::Data& data, 
					   const Data::Fields& fields, 
					   const AmanziGeometry::GeometricModelPtr& gm)
{
  const int num_element_blocks = data.element_blocks ();
  const int num_side_sets      = data.side_sets ();
  const int num_node_sets      = data.node_sets ();

  const int space_dimension = data.parameters ().dimensions ();

  // Something to put elements in

  elements_part_ = &meta_data_->declare_part ("Elements", element_rank_);

  // Convert element blocks, node and side sets into Parts:

  for (int block = 0; block < num_element_blocks; ++block) {
    element_blocks_.push_back (add_element_block_ (data.element_block (block)));
  }

  // Something to put nodes in
  nodes_part_ = &meta_data_->declare_part ("Nodes", node_rank_);

  // Add a faces part.
  faces_part_ = &meta_data_->declare_part ("Element sides", face_rank_);

  // Add a field to faces that represents the face "owner"
    
  face_owner_ = &( meta_data_->declare_field< Id_field_type >("FaceOwner") );
  stk::mesh::put_field(*face_owner_, face_rank_, *faces_part_);

  for (int side_set = 0; side_set < num_side_sets; ++side_set)
    side_sets_.push_back (add_side_set_ (data.side_set (side_set)));

  for (int node_set = 0; node_set < num_node_sets; ++node_set)
    node_sets_.push_back (add_node_set_ (data.node_set (node_set)));

  // Get the universal part. There's only one everything.
  stk::mesh::Part &universal_part (meta_data_->universal_part ());

  put_coordinate_field_ (universal_part, space_dimension);

  // Declare and Put fields
  for (Data::Fields::const_iterator field = fields.begin ();
       field != fields.end (); 
       ++field)
  {
    put_field_ (*field, universal_part, space_dimension);
  }


  // Add some additional parts as place holders for meshets that will
  // be created based on the region specifications in the geometric model

  if (gm)
    init_extra_parts_from_gm(gm);

  meta_data_->commit ();

}

/** 
 * add contents of the blocks and sets to their respctive parts
 * 
 * @param data 
 * @param cellmap 
 * @param vertmap 
 * @param fields 
 */
void Mesh_STK_factory::build_bulk_data_ (const Data::Data& data, 
					 const Epetra_Map& cellmap, 
					 const Epetra_Map& vertmap,
					 const Data::Fields& fields,
					 const AmanziGeometry::GeometricModelPtr& gm)
{
  const int space_dimension = data.parameters ().dimensions ();

  // elements first, nodes will be declared and will overlap on processors

  bulk_data_->modification_begin();

  unsigned int localcidx(0);
  for (unsigned int block = 0; block < element_blocks_.size (); ++block) {
    add_elements_to_part_ (data.element_block (block), *element_blocks_ [block], 
                           cellmap, vertmap, localcidx);
  }
    
  bulk_data_->modification_end();

  // Now node ownership is established, so we can start doing things
  // with nodes; bad things happen if you try to do these things
  // before ownership is established

  bulk_data_->modification_begin();

  for (int set = 0; set < node_sets_.size (); ++set)
    add_nodes_to_part_ (data.node_set (set), *node_sets_ [set], vertmap);
  add_coordinates_ (data.coordinates (), vertmap);

  bulk_data_->modification_end();

  // Now, generate unique faces.  I have not found a way to renumber
  // entities after they are declared.  So, we need to know the
  // global face id in advance.  That means they need to be counted
  // and starting global indexes distributed.
    
  int nface_local(count_local_faces_());
  std::vector<int> nface(communicator_.NumProc(), 0);
  communicator_.GatherAll(&nface_local, &nface[0], 1);

  std::vector<int> global_fidx(nface.size() + 1);
  int nface_global(0);
  global_fidx[0] = 1;
  for (unsigned int i = 1; i < global_fidx.size(); i++) {
    global_fidx[i] = global_fidx[i-1] + nface[i-1];
    nface_global += nface[i-1];
  }

  // std::cerr << "Process " << communicator_.MyPID() 
  //           << " owns " << nface_local 
  //           << " of " << nface_global 
  //           << " faces, index " << global_fidx[communicator_.MyPID()]
  //           << " to "  << global_fidx[communicator_.MyPID()+1] - 1
  //           << std::endl;
  ASSERT((global_fidx.back()-1) == nface_global);

  bulk_data_->modification_begin();
  generate_local_faces_(global_fidx[communicator_.MyPID()], false);
  bulk_data_->modification_end();

  // Put side set faces in the correct parts
    
  check_face_ownership_();
  generate_cell_face_relations();
  bulk_data_->modification_begin();
  for (int set = 0; set < side_sets_.size (); ++set) {
    add_sides_to_part_ (data.side_set (set), *side_sets_ [set], cellmap);
  }
  bulk_data_->modification_end();

  // Fill in the elements/faces/nodes for the extra parts that are created
  // based on the region specifications in the geometric model

  if (gm)
    fill_extra_parts_from_gm(gm);
}


void Mesh_STK_factory::put_field_ (const Data::Field& field_data, 
                               stk::mesh::Part& part, 
                               unsigned int space_dimension)
{
  const unsigned int location = entity_map_->kind_to_rank (field_data.location ());

  if (field_data.type () == Data::SCALAR)
  {
    Scalar_field_type& field (meta_data_->
                              declare_field<Scalar_field_type>(field_data.name ()));
    stk::mesh::put_field (field, location, part);
  }

  if (field_data.type () == Data::VECTOR)
  {
    Vector_field_type& field (meta_data_->
                              declare_field<Vector_field_type>(field_data.name ()));
    stk::mesh::put_field (field, location, part, space_dimension);
  }

}

void Mesh_STK_factory::put_coordinate_field_ (stk::mesh::Part& part, unsigned int space_dimension)
{
  coordinate_field_ = & meta_data_->declare_field<Vector_field_type>("coordinates");
  stk::mesh::put_field (*coordinate_field_, stk::mesh::Node, part, space_dimension);
}

/** 
 * This routine puts the coordinate fields on @em local nodes.  
 *
 * Caller is responsible for putting the ::bulk_data_ into
 * modification state.
 * 
 * The specified @c coordinate_data is for the nodes involved in cells
 * owned by the local processor.  
 * 
 * @param coordinate_data coordinates of nodes originally declared on this processor
 * @param vertmap local to global node index map for nodes originally declared on this processor
 */
void Mesh_STK_factory::add_coordinates_ (const Data::Coordinates<double>& coordinate_data,
                                     const Epetra_Map& vertmap)
{

  // Select all the local nodes this process knows about
  stk::mesh::Selector owned(meta_data_->universal_part ());
  Entity_vector local_nodes;
  stk::mesh::get_selected_entities (owned,
                                    bulk_data_->buckets (stk::mesh::Node), 
                                    local_nodes);

  // Loop over the nodes, if the node is used by the cells owned by
  // this process, set the coordinate
  int node_coordinate_index = 0;
  for (Entity_vector::const_iterator node_it = local_nodes.begin ();
       node_it != local_nodes.end (); ++node_it)
  {
    int global_vidx((*node_it)->identifier());
    if (vertmap.MyGID(global_vidx)) {
      // only set the coordinates for the nodes this process
      // knows about
      int local_vidx(vertmap.LID(global_vidx));
      double * coordinate_field_data = 
          stk::mesh::field_data (*coordinate_field_, **node_it);
      coordinate_data (local_vidx, coordinate_field_data);
      // std::cerr << "add_coordinates: node " << global_vidx << " (" 
      //           << local_vidx << " local): " 
      //           << coordinate_field_data[0] << ", "
      //           << coordinate_field_data[1] << ", "
      //           << coordinate_field_data[2] << ", "
      //           << std::endl;
    }
  }
}


  // Declare a part from an element block

stk::mesh::Part* Mesh_STK_factory::add_element_block_ (const Data::Element_block& block)
{

  std::ostringstream name;
  if (block.name ().size () > 0)
    name << block.name ();
  else
    name << "element block " << block.id ();

  stk::mesh::Part &new_part (meta_data_->declare_part (name.str (), element_rank_));

  meta_data_->declare_part_subset(*elements_part_, new_part);

  Cell_type mdtype(block.element_type ());
  switch (mdtype) {
    case HEX:
    case TET:
    case PYRAMID:
    case PRISM:
      // this type is OK
      break;
    default:
      std::string msg = 
          boost::str(boost::format("%s has unsupported cell type %s") %
                     name.str() % Data::type_to_name(block.element_type()));
      Exceptions::amanzi_throw( STK::Error (msg.c_str()) );
      break;
  }

  const shards::CellTopology topo(get_topology_data(mdtype));
      
  stk::mesh::set_cell_topology (new_part, topo.getTopology ());

  add_set_part_relation_ (block.id (), new_part);

  return &new_part;
}


  // Declare a part from just a name and id (we don't yet know the elements that will go into it)

  stk::mesh::Part* Mesh_STK_factory::add_element_block_ (const std::string name, const int id)
{

  stk::mesh::Part &new_part (meta_data_->declare_part (name, element_rank_));

  meta_data_->declare_part_subset(*elements_part_, new_part);

  add_set_part_relation_ (id, new_part);

  return &new_part;
}


  // Declare a part from a side set data structure

stk::mesh::Part* Mesh_STK_factory::add_side_set_ (const Data::Side_set& set)
{
  std::ostringstream name;
  if (set.name ().size () > 0)
    name << set.name ();
  else
    name << "side set " << set.id ();
  stk::mesh::Part &new_part (meta_data_->declare_part (name.str (), face_rank_));
  meta_data_->declare_part_subset(*faces_part_, new_part);

  add_set_part_relation_ (set.id (), new_part);

  return &new_part;
}

  // Declare a part from a side set name and id

  stk::mesh::Part* Mesh_STK_factory::add_side_set_ (const std::string name, const int id)
{
  stk::mesh::Part &new_part (meta_data_->declare_part (name, face_rank_));
  meta_data_->declare_part_subset(*faces_part_, new_part);

  add_set_part_relation_ (id, new_part);

  return &new_part;
}


  // Declare a part from a node set data structure

stk::mesh::Part* Mesh_STK_factory::add_node_set_ (const Data::Node_set& set)
{
  std::ostringstream name;
  if (set.name ().size () > 0)
    name << set.name ();
  else
    name << "node set " << set.id ();
  stk::mesh::Part &new_part (meta_data_->declare_part (name.str (), stk::mesh::Node));
  meta_data_->declare_part_subset(*nodes_part_, new_part);

  add_set_part_relation_ (set.id (), new_part);

  return &new_part;
}

  stk::mesh::Part* Mesh_STK_factory::add_node_set_ (const std::string name, const int id)
{
  stk::mesh::Part &new_part (meta_data_->declare_part (name, stk::mesh::Node));
  meta_data_->declare_part_subset(*nodes_part_, new_part);

  add_set_part_relation_ (id, new_part);

  return &new_part;
}

void Mesh_STK_factory::add_set_part_relation_ (unsigned int set_id, stk::mesh::Part& part)
{
  const unsigned int part_id = part.mesh_meta_data_ordinal ();
  const stk::mesh::EntityRank rank = part.primary_entity_rank ();
  const Rank_and_id rank_set_id = std::make_pair (rank, set_id);

  ASSERT (set_to_part_.find (rank_set_id) == set_to_part_.end ());

  set_to_part_ [rank_set_id]  = &part;

}


/** 
 * This routine adds elements to the specified part.  Each element is
 * expected to have the same topology and this topology must match
 * that assigned to the part.  In the process, nodes will be declared,
 * but node ownership will not be assigned.
 * 
 * @param block 
 * @param part 
 * @param cmap 
 * @param vmap 
 * @param localidx0 current local element index (updated)
 */
void Mesh_STK_factory::add_elements_to_part_ (const Data::Element_block& block, stk::mesh::Part &part,
                                          const Epetra_Map& cmap, const Epetra_Map& vmap, unsigned int& localidx0)
{
  // Add connectivity information via stk::mesh::declare_element
  std::vector<int> storage (block.nodes_per_element ());
  std::vector<int> global_vidx(block.nodes_per_element ());

  for (int bidx = 0; bidx < block.num_elements (); ++bidx, ++localidx0)
  {
    block.connectivity (bidx, storage.begin ());

    for (unsigned int i = 0; i < block.nodes_per_element (); i++) {
      global_vidx[i] = vmap.GID(storage[i]);
    }

    int global_cidx(cmap.GID(localidx0));

    try {
      stk::mesh::Entity &element = 
          stk::mesh::declare_element (*bulk_data_, part, global_cidx, &global_vidx[0]);
    } catch (const std::exception& e) {
      std::stringstream msg;
      std::cerr << e.what() << std::endl;
      msg << "cell error: local: " << localidx0 << ": ";
      std::copy(storage.begin(), storage.end(), std::ostream_iterator<int>(msg, ", "));
      msg << "global: " << global_cidx << ": ";
      std::copy(global_vidx.begin(), global_vidx.end(), std::ostream_iterator<int>(msg, ", "));
      std::cerr << msg.str() << std::endl;

      Exceptions::amanzi_throw( STK::Error (msg.str().c_str()) );
    }

  }
}

Entity_vector
Mesh_STK_factory::get_element_side_nodes_(const stk::mesh::Entity& element, 
                                      const unsigned int& s)
{
  ASSERT(element.entity_rank() == element_rank_);

  const CellTopologyData* topo = stk::mesh::get_cell_topology (element);
  ASSERT(topo != NULL);

  // get the cell nodes on this side

  const CellTopologyData *side_topo = (topo->side[s].topology);
  const unsigned * const side_node_map = topo->side[s].node;

  stk::mesh::PairIterRelation rel = element.relations( node_rank_ );
  Entity_vector snodes;
  for ( unsigned i = 0 ; i < side_topo->node_count ; ++i ) {
    snodes.push_back(rel[ side_node_map[i] ].entity());
  }
  return snodes;
}



/** 
 * This routine finds (indirectly through related nodes) the face
 * declared for the specified element side.  
 * 
 * @param element The element
 * @param s The side index (0-based)
 * 
 * @return face entity, or NULL if non-existent
 */
stk::mesh::Entity *
Mesh_STK_factory::get_element_side_face_(const stk::mesh::Entity& element, const unsigned int& s)
{
  Entity_vector snodes(get_element_side_nodes_(element, s));

  // get the face the nodes relate to, if any

  std::vector< stk::mesh::Entity * > rfaces;
  get_entities_through_relations(snodes, face_rank_, rfaces);
  
  ASSERT(rfaces.size() < 2);

  stk::mesh::Entity *result = NULL;
  if (!rfaces.empty()) result = rfaces.front();

  return result;
}
    



/** 
 * This routine takes care of declaring a (local) face, include
 * declaring the entity and making necessary relations to related
 * nodes.  Cell relations are declared in
 * ::generate_cell_face_relations().
 * 
 * @param nodes nodes involved in the face
 * @param index global index to assign to the face
 * 
 * @return the declared face entity
 */
const stk::mesh::Entity&
Mesh_STK_factory::declare_face_(stk::mesh::EntityVector& nodes, 
                            const unsigned int& index, 
                            const unsigned int& owner_index)
{
  stk::mesh::PartVector p;
  p.push_back(faces_part_);

  stk::mesh::Entity& face = 
      bulk_data_->declare_entity(face_rank_, index, p);

  stk::mesh::FieldTraits<Id_field_type>::data_type *owner = 
      stk::mesh::field_data(*face_owner_, face);
  *owner = owner_index;

  unsigned int r(0);
  for (stk::mesh::EntityVector::iterator n = nodes.begin(); n != nodes.end(); n++) {
    bulk_data_->declare_relation(face, **n, r++);
  }

  return face;
}

/** 
 * Local
 * 
 * This routine goes through the locally owned cells and generates
 * cell faces.  The result are faces that are uniquely owned and
 * uniquely defined.  For example, only one face relating two
 * neighboring cells will be declared.
 *
 * Some conventions:
 *
 *  - A face will be declared by the process that owns the related
 *    cell(s).  N.B.: this can be changed by stk::mesh
 * 
 *  - If the processes owning cells related to a given face are
 *    different, the face is declared by the lower ranked process.
 *
 * This routine just declares the needed faces.  It would be nice if
 * stk::mesh would maintain ownership as the declaration, but, of
 * course it doesn't.  In order to deal with side sets, ownership
 * needs to be established/changed (done in ::check_face_ownership_()).
 *
 * caller is responsible for calling
 * stk::mesh::BulkData::modification_begin() (when @c justcount is @c
 * false ).
 *
 * This would be a lot easier if we could just declare the faces with
 * some bogus global indexes and then renumber them, but I can't
 * figure out how to do that.  So, the needed faces have to be counted
 * first to establish global index ranges for each processor.
 *
 * Face to node relations are done here when a face is declared.  Cell
 * to face relations are done afterward when declared faces are
 * available to all processes (see ::generate_cell_face_relations())
 * 
 * @param faceidx0 starting global face index
 *
 * @param justcount if true, do not declare any faces, just count how
 * many there would be.
 * 
 * @return number of faces declared/counted
 */
int 
Mesh_STK_factory::generate_local_faces_(const int& faceidx0, const bool& justcount)
{      
  const unsigned int me = communicator_.MyPID();
  int faceidx(faceidx0);
  stk::mesh::Selector owned(meta_data_->locally_owned_part());

  int nlocal = 
      stk::mesh::count_selected_entities(owned, bulk_data_->buckets(element_rank_));

  stk::mesh::PartVector parts(meta_data_->get_parts());
  for (stk::mesh::PartVector::iterator p = parts.begin(); 
       p != parts.end(); p++) {
    const CellTopologyData* topo = stk::mesh::get_cell_topology (**p);
        
    if (topo == NULL) continue;

    stk::mesh::Selector s(*(*p));
    s &= owned;

    std::vector< stk::mesh::Entity * > cells;
    stk::mesh::get_selected_entities(s, 
                                     bulk_data_->buckets(element_rank_), 
                                     cells);
        
    int localidx(0);
    for (std::vector< stk::mesh::Entity * >::iterator c = cells.begin();
         c != cells.end(); c++, localidx++) {

      for (unsigned int s = 0; s < topo->side_count; s++) {

        // see if the face already exists (on the local
        // processor); if so, we don't need to bother with

        stk::mesh::Entity *rface(get_element_side_face_(**c, s));
        if (rface != NULL) { 
          // std::cerr << 
          //   boost::str(boost::format("%d: error: cell %d, side %d: will not make duplicate face") %
          //              me % (*c)->identifier() % s) << std::endl;
          continue;
        }
                                         
        // get the cell nodes on this side

        Entity_vector rcells;
        Entity_vector snodes(get_element_side_nodes_(**c, s));

        // get whatever cells the nodes on this side are
        // related to (*c should be one of them, *and* there
        // should be at most one other cell)
            
        get_entities_through_relations(snodes, 
                                       element_rank_, 
                                       rcells);

        if (rcells.size() < 1) {
                    
          std::string msg = 
              boost::str(boost::format("error: cell %d (owner = %d): face %d (%d nodes) doesn't even relate to cell") %
                         (*c)->identifier() % me % s % snodes.size());
          // std::cerr << msg << std::endl;
          Exceptions::amanzi_throw( STK::Error (msg.c_str()) );
        } else if (rcells.size() > 2) {
          std::string msg = 
              boost::str(boost::format("error: cell %d (owner = %d): face %d relates to too many cells (%d)") %
                         (*c)->identifier() % me % s % rcells.size());
          Exceptions::amanzi_throw( STK::Error (msg.c_str()) );
        }

        stk::mesh::Entity *rcell(NULL);

        for ( std::vector< stk::mesh::Entity * >::iterator r = rcells.begin();
              r != rcells.end(); r++) {
          if ((*r != *c)) {
            rcell = *r;
          }
        }


        // if there is no cell on the other side (rcell ==
        // NULL) of the face, just declare/count the face

        // if the related cell is owned by another
        // process, only make the face if this processor's
        // rank is less; 

        // if the related cell is owned by this process, only make it
        // for the cell with the lower id
                
        if (rcell != NULL) {
          if (rcell->owner_rank() < me) {
            continue;
          } else if (rcell->owner_rank() == me) {
            if (rcell->identifier() < (*c)->identifier()) {
              continue;
            }
          }
        }


        if (!justcount) {
          const stk::mesh::Entity& face = 
            declare_face_(snodes, faceidx, (*c)->identifier());
        }
        faceidx++;
                
        // std::cerr << "Cell " << (*c)->identifier() 
        //           << ": side " << s 
        //           << ": face idx = " << faceidx
        //           << std::endl;

      } // side loop

    } // entity (element) loop
 
  } // part loop
    
  return faceidx - faceidx0;
}

// -------------------------------------------------------------
// Mesh_STK_factory::generate_cell_face_relations
// -------------------------------------------------------------
/** 
 * Collective
 *
 * This routine generates cell-to-face relations.  It must be called
 * after ::generate_local_faces_().  It needs to be a separate
 * transaction so that all processes are aware of declared faces and
 * cells.
 * 
 */
void
Mesh_STK_factory::generate_cell_face_relations(void)
{
  const unsigned int me = communicator_.MyPID();
  stk::mesh::Selector owned(meta_data_->locally_owned_part());
  stk::mesh::PartVector parts(meta_data_->get_parts());

  bulk_data_->modification_begin();

  for (stk::mesh::PartVector::iterator p = parts.begin(); 
       p != parts.end(); p++) {

    const CellTopologyData* topo = stk::mesh::get_cell_topology (**p);
        
    if (topo == NULL) continue;

    stk::mesh::Selector s(*(*p));
    s &= owned;

    std::vector< stk::mesh::Entity * > cells;
    stk::mesh::get_selected_entities(s, 
                                     bulk_data_->buckets(element_rank_), 
                                     cells);
        
    int localidx(0);
    for (std::vector< stk::mesh::Entity * >::iterator c = cells.begin();
         c != cells.end(); c++, localidx++) {

      for (unsigned int s = 0; s < topo->side_count; s++) {

        // see if the face already exists (on the local
        // processor); if so, see if this cell is already
        // related to it (it shouldn't be)

        stk::mesh::Entity *rface(get_element_side_face_(**c, s));
        bulk_data_->declare_relation(*(*c), *rface, s);
      }
    }
  }

  bulk_data_->modification_end();
}

/** 
 * Collective
 *
 * It appears that stk::mesh requires that the local process @em own
 * an entity in order to call stk::mesh::change_entity_parts() for it.
 * Sometimes, but not always, stk::mesh assigns ownership of faces
 * different than they were declared.  This probably has something to
 * do with the ownership of the nodes involved.  
 *
 * More testing/knowledge is needed to be certain.  
 *
 * In order to deal with side sets, as stored in Data::Side_set,
 * the face must be owned by the same process as the cell in which it
 * is involved.  This routine makes sure that's the case.
 *
 * The convention for faces that connect cells is for the process with
 * lowest rank to own the face.
 * 
 */
void
Mesh_STK_factory::check_face_ownership_(void)
{
  stk::mesh::check_face_ownership(*bulk_data_);
}

// -------------------------------------------------------------
// Mesh_STK_factory::add_sides_to_part_
// -------------------------------------------------------------

void Mesh_STK_factory::add_sides_to_part_ (const Data::Side_set& side_set, 
                                           stk::mesh::Part &part,
                                           const Epetra_Map& cmap)
{
  const unsigned int me(communicator_.MyPID());

  // Side set consists of elements (local index) and a local face
  // number. We need to convert these to the unique face indices.

  stk::mesh::PartVector parts_to_add;
  parts_to_add.push_back (&part);

  const int num_sides  = side_set.num_sides ();
  const std::vector<int>& element_list = side_set.element_list ();
  const std::vector<int>& side_list    = side_set.side_list ();
  ASSERT (element_list.size () == num_sides);
  ASSERT (side_list.size ()    == num_sides);

  for (int index=0; index < num_sides; ++index)
  {
    int local_idx(element_list [index]);
    ASSERT(local_idx < cmap.NumMyElements());

    const int local_side = side_list [index];
    int global_idx(cmap.GID(element_list [index]));

    // Look up the element from the Id.
    stk::mesh::Entity *element = bulk_data_->get_entity (element_rank_, global_idx);

    stk::mesh::Entity *face = get_element_side_face_(*element, local_side);

    if (face == NULL) {
      std::string msg = 
          boost::str(boost::format("side set %d: face not declared for side %d of cell %d (local = %d") %
                     side_set.id() % local_side % global_idx % local_idx);
      Exceptions::amanzi_throw( STK::Error (msg.c_str()) );
    }

    // std::cerr << 
    //   boost::str(boost::format("%02d: Side set %2d: %5d (%5d) %2d (%5d)") % 
    //              me % side_set.id() % global_idx % local_idx % local_side % face->identifier())
    //           << std::endl;

    if (face->owner_rank() != me) {
      std::string msg = 
          boost::str(boost::format("Process %d: side set %d: I own cell %d, but not its side face %d") %
                     me % side_set.id() % element->identifier() % face->identifier());
      std::cerr << msg << std::endl;
    } else {
      bulk_data_->change_entity_parts(*face, parts_to_add);
    }
  }
}


/** 
 * This routine puts the nodes in the node set in the specified part.
 * Ownership of nodes should be established before this is
 * called. Only those nodes owned by this process are put in the part.
 * See comments for ::add_coordinates_() for an explanation of
 * indexing.
 * 
 * @param node_set 
 * @param part 
 * @param vmap 
 */
void Mesh_STK_factory::add_nodes_to_part_ (const Data::Node_set& node_set, 
                                       stk::mesh::Part &part,
                                       const Epetra_Map& vmap)
{
  const unsigned int me(communicator_.MyPID());

  stk::mesh::PartVector parts_to_add;
  parts_to_add.push_back (&part);

  const int num_nodes = node_set.num_nodes ();
  const std::vector<int>& node_list = node_set.node_list ();
  ASSERT (node_list.size () == num_nodes);

  for (std::vector<int>::const_iterator it = node_list.begin ();
       it != node_list.end ();
       ++it)
  {
    int local_vidx(*it);
    int global_vidx(vmap.GID(local_vidx));
    stk::mesh::Entity *node = bulk_data_->get_entity (node_rank_, global_vidx);
    if (node->owner_rank() == me) {
      bulk_data_->change_entity_parts(*node, parts_to_add);
    }
  }

}





void Mesh_STK_factory::init_extra_parts_from_gm(const AmanziGeometry::GeometricModelPtr& gm) 
{

  if (gm == NULL) return;


  int space_dim = gm->dimension();


  int ngr = gm->Num_Regions();
  for (int i = 0; i < ngr; i++)
    {
      AmanziGeometry::RegionPtr greg = gm->Region_i(i);
      
      switch (greg->type())
        {
        case AmanziGeometry::BOX: { 

          if (((AmanziGeometry::BoxRegionPtr) greg)->is_degenerate()) {            
            add_side_set_(greg->name(), greg->id());
          }
          else {
            add_element_block_(greg->name(), greg->id());
          }

          break;
        }
        case AmanziGeometry::PLANE: {

          // Assumption is that user wants to extract side sets
          
          add_side_set_(greg->name(), greg->id());

          break;
        }
        case AmanziGeometry::LABELEDSET: {

          // This should already be in there

          break;
          }
        case AmanziGeometry::SURFACE: {
          // Not implemented


          break;
        }
        case AmanziGeometry::POINT: {

          // Cell set based on points

          add_element_block_(greg->name(), greg->id());

          break;
        }
        case AmanziGeometry::COLORFUNCTION: {

	  // Cell set based on points

          add_element_block_(greg->name(), greg->id());

          break;
        }
        default:
          throw std::exception();
        }
      
    }

}




void Mesh_STK_factory::fill_extra_parts_from_gm(const AmanziGeometry::GeometricModelPtr& gm) 
{
  if (gm == NULL) return;


  bulk_data_->modification_begin();

  int space_dim = gm->dimension();


  int ngr = gm->Num_Regions();
  for (int i = 0; i < ngr; i++)
    {
      AmanziGeometry::RegionPtr greg = gm->Region_i(i);
      stk::mesh::Part *part;
      stk::mesh::PartVector parts_to_add;

      switch (greg->type())
        {
        case AmanziGeometry::BOX: { 

          // Find the part with this name

	  part = meta_data_->get_part (greg->name());
          ASSERT (part);
	  parts_to_add.push_back(part);

          if (((AmanziGeometry::BoxRegionPtr) greg)->is_degenerate()) {

            ASSERT (part->primary_entity_rank () == face_rank_);

            // Assumption is that user wants to extract cell sets
          
            stk::mesh::Selector owned(meta_data_->locally_owned_part());
            stk::mesh::EntityVector faces;
            stk::mesh::get_selected_entities(owned, bulk_data_->buckets(stk::mesh::Face), faces);

            stk::mesh::EntityVector::iterator f;
            for (f = faces.begin(); f != faces.end(); f++) {

              stk::mesh::PairIterRelation nodes = (*f)->relations (node_rank_);

              double cen[3]={0.0,0.0,0.0};
	      int nfn = 0;
              for (stk::mesh::PairIterRelation::iterator it = nodes.begin (); 
                   it != nodes.end (); ++it)
                {
                  double *xyz = stk::mesh::field_data(*coordinate_field_, *(it->entity()));
                  for (int k = 0; k < space_dim; k++)
                    cen[k] += xyz[k];
		  nfn++;
                }
              for (int k = 0; k < space_dim; k++)
                cen[k] /= nfn;

              AmanziGeometry::Point pcen(space_dim);
              pcen.set(cen);
              
              if (greg->inside(pcen))  // If face center is inside region
                {
                  bulk_data_->change_entity_parts(*(*f),parts_to_add);
                }
            }

          }

          else {

            // Assumption is that user wants to extract cell sets

            ASSERT (part->primary_entity_rank () == element_rank_);

          
            stk::mesh::Selector owned(meta_data_->locally_owned_part());
            stk::mesh::EntityVector cells;
            stk::mesh::get_selected_entities(owned, bulk_data_->buckets(stk::mesh::Element), cells);

            stk::mesh::EntityVector::iterator c;
            for (c = cells.begin(); c != cells.end(); c++) {

              stk::mesh::PairIterRelation nodes = (*c)->relations (node_rank_);

              double cen[3]={0.0,0.0,0.0};
	      int nen = 0;
              for (stk::mesh::PairIterRelation::iterator it = nodes.begin (); 
                   it != nodes.end (); ++it)
                {
                  double *xyz = stk::mesh::field_data(*coordinate_field_, *(it->entity()));
                  for (int k = 0; k < space_dim; k++)
                    cen[k] += xyz[k];
		  nen++;
                }
              for (int k = 0; k < space_dim; k++)
                cen[k] /= nen;

              AmanziGeometry::Point pcen(space_dim);
              pcen.set(cen);
              
              if (greg->inside(pcen))  // If center is inside region
                {
                  bulk_data_->change_entity_parts(*(*c),parts_to_add);
                }
            }

          }

          break;
        }
        case AmanziGeometry::PLANE: {

	  part = meta_data_->get_part (greg->name());
          ASSERT (part);
	  parts_to_add.push_back(part);
	  
          // Assumption is that user wants to extract side sets
	  
	  ASSERT (part->primary_entity_rank () == face_rank_);

	  stk::mesh::Selector owned(meta_data_->locally_owned_part());
	  stk::mesh::EntityVector faces;
	  stk::mesh::get_selected_entities(owned, bulk_data_->buckets(stk::mesh::Face), faces);
	  
	  stk::mesh::EntityVector::iterator f;
	  for (f = faces.begin(); f != faces.end(); f++) {
	    
	    stk::mesh::PairIterRelation nodes = (*f)->relations (node_rank_);
	    
	    bool on_plane = true;
	    for (stk::mesh::PairIterRelation::iterator it = nodes.begin (); 
		 it != nodes.end (); ++it)
	      {
		double *xyz = stk::mesh::field_data(*coordinate_field_, *(it->entity()));
		AmanziGeometry::Point pxyz(space_dim);
		pxyz.set(xyz);
		
		if (!greg->inside(pxyz)) {
		  on_plane = false;
		  break;
		}
	      }
	    
	    if (on_plane)  // If face center is inside region
	      {
		bulk_data_->change_entity_parts(*(*f),parts_to_add);
	      }
	  }

          break;
        }
        case AmanziGeometry::LABELEDSET: {

          // This should already be in there

          break;
	}
        case AmanziGeometry::SURFACE: {
          // Not implemented


          break;
        }
        case AmanziGeometry::POINT: {

	  part = meta_data_->get_part (greg->name());
          ASSERT (part);
	  parts_to_add.push_back(part);
	  
	  
	  stk::mesh::Selector owned(meta_data_->locally_owned_part());
	  stk::mesh::EntityVector nodes;
	  stk::mesh::get_selected_entities(owned, bulk_data_->buckets(stk::mesh::Node), nodes);
	 

          AmanziGeometry::Point rgnpnt(((AmanziGeometry::PointRegionPtr)greg)->point());

          double mindist2 = 1.0e+16;
          stk::mesh::Entity *minnode = NULL;
 
	  stk::mesh::EntityVector::iterator n;
	  for (n = nodes.begin(); n != nodes.end(); n++) {
	    
            double *xyz = stk::mesh::field_data(*coordinate_field_, *(*n));
            AmanziGeometry::Point vpnt(space_dim);
            vpnt.set(xyz);

            double dist2 = (vpnt-rgnpnt)*(vpnt-rgnpnt);

            if (dist2 < mindist2) {              
              mindist2 = dist2;
              minnode = *n;
              if (mindist2 <= 1.0e-32)
                break;
            }            
		
	  }

          stk::mesh::PairIterRelation cells = minnode->relations (element_rank_);
	    
          for (stk::mesh::PairIterRelation::iterator it = cells.begin(); 
               it != cells.end(); it++) {

            stk::mesh::Entity *element = it->entity();


            // build up a description of the element in terms of its nodes ...

            std::vector<AmanziGeometry::Point> ccoord;

            stk::mesh::PairIterRelation elnodes = element->relations(node_rank_);
            for (stk::mesh::PairIterRelation::iterator jt = elnodes.begin();
                 jt != elnodes.end(); jt++) {

              stk::mesh::Entity *elnode = jt->entity();

              double *xyz = stk::mesh::field_data(*coordinate_field_, *elnode);
              AmanziGeometry::Point vpnt(space_dim);
              vpnt.set(xyz);
              ccoord.push_back(vpnt);
              
            }

            
            // ... and its faces (while ensuring that they are pointing out of
            // the element)

            std::vector<AmanziGeometry::Point> fcoord;
            std::vector<unsigned int> nfnodes;

            stk::mesh::PairIterRelation elfaces = element->relations(face_rank_);

            int nf = elfaces.size();

            for (stk::mesh::PairIterRelation::iterator jt = elfaces.begin (); 
                 jt != elfaces.end (); ++jt)
              {
                stk::mesh::Entity *elface = jt->entity();
                
                stk::mesh::FieldTraits<Id_field_type>::data_type *owner = 
                  stk::mesh::field_data<Id_field_type>(*face_owner_, *elface);
                
                int dir = (*owner == element->identifier()) ? 1 : -1;
                
                
                stk::mesh::PairIterRelation fnodes = elface->relations(node_rank_);                
                nfnodes.push_back(fnodes.size());
                
                std::vector<AmanziGeometry::Point> thisfcoord;
                for (stk::mesh::PairIterRelation::iterator kt = fnodes.begin();
                     kt != fnodes.end(); ++kt)
                  {
                    stk::mesh::Entity *fnode = kt->entity();
                    
                    double *xyz = stk::mesh::field_data(*coordinate_field_, *fnode);
                    AmanziGeometry::Point vpnt(space_dim);
                    vpnt.set(xyz);
                    thisfcoord.push_back(vpnt);
                  }       
                
                if (dir == 1) {
                  for (std::vector<AmanziGeometry::Point>::iterator kt = thisfcoord.begin(); kt != thisfcoord.end(); ++kt)
                    {
                      fcoord.push_back(*kt);
                    }
                }
                else {
                  for (std::vector<AmanziGeometry::Point>::reverse_iterator kt = thisfcoord.rbegin(); kt != thisfcoord.rend(); ++kt)
                    {
                      fcoord.push_back(*kt);
                    }
                }
              }
                        

            // Test if the element polyhedron contains the rgn pnt

            if (AmanziGeometry::point_in_polyhed(rgnpnt,ccoord,nf,nfnodes,
                                                       fcoord)) {
              bulk_data_->change_entity_parts(*element,parts_to_add);
            }

          }

          break;
        }
        case AmanziGeometry::COLORFUNCTION: { 

          // Find the part with this name
          
	  part = meta_data_->get_part (greg->name());
          ASSERT (part);
	  parts_to_add.push_back(part);
          
          
          // Assumption is that user wants to extract cell sets
          
          ASSERT (part->primary_entity_rank () == element_rank_);
          
          
          stk::mesh::Selector owned(meta_data_->locally_owned_part());
          stk::mesh::EntityVector cells;
          stk::mesh::get_selected_entities(owned, bulk_data_->buckets(stk::mesh::Element), cells);
          
          stk::mesh::EntityVector::iterator c;
          for (c = cells.begin(); c != cells.end(); c++) {
            
            stk::mesh::PairIterRelation nodes = (*c)->relations (node_rank_);
            
            double cen[3]={0.0,0.0,0.0};
            int nen = 0;
            for (stk::mesh::PairIterRelation::iterator it = nodes.begin (); 
                 it != nodes.end (); ++it)
              {
                double *xyz = stk::mesh::field_data(*coordinate_field_, *(it->entity()));
                for (int k = 0; k < space_dim; k++)
                  cen[k] += xyz[k];
                nen++;
              }
            for (int k = 0; k < space_dim; k++)
              cen[k] /= nen;
            
            AmanziGeometry::Point pcen(space_dim);
            pcen.set(cen);
            
            if (greg->inside(pcen))  // If center is inside region
              {
                bulk_data_->change_entity_parts(*(*c),parts_to_add);
              }
          }
          
          break;
        }
        default:
          throw std::exception();
        }
      
    }

}



} // close namespace STK 
} // close namespace Mesh 
} // close namespace Amanzi 


