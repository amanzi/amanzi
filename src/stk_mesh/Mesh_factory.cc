#include "Mesh_factory.hh"
#include "Mesh.hh"

#include "dbc.hh"

#include <iostream>
#include <algorithm>
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


// STK_mesh
#include "Element_field_types.hh"
#include "Element_category.hh"
#include "Cell_topology.hh"

// Trilinos
#include <Shards_CellTopology.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologicalMetaData.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetEntities.hpp>

namespace STK_mesh
{

Mesh_factory::Mesh_factory (const stk::ParallelMachine& comm, int bucket_size) 
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
Mesh* Mesh_factory::build_mesh (const Mesh_data::Data& data, 
                                const Mesh_data::Fields& fields)
{
    ASSERT(communicator_.NumProc() == 1);

    int ncell(data.parameters().num_elements_);
    Epetra_Map cmap(ncell, 1, communicator_);

    int nvert(data.parameters().num_nodes_);
    Epetra_Map vmap(nvert, 1, communicator_);

    return build_mesh(data, cmap, vmap, fields);
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
Mesh* Mesh_factory::build_mesh (const Mesh_data::Data& data, 
                                const Epetra_Map& cellmap,
                                const Epetra_Map& vertmap,
                                const Mesh_data::Fields& fields)
{

    // Update construction variables for the new mesh.
    const int space_dimension = data.parameters ().dimensions ();

    ASSERT (Mesh::valid_dimension (space_dimension));

    entity_map_ = new Entity_map (space_dimension);
    meta_data_  = new stk::mesh::MetaData (stk::mesh::fem_entity_rank_names ());
    bulk_data_  = new stk::mesh::BulkData (*meta_data_, parallel_machine_, bucket_size_);

    node_rank_    = entity_map_->kind_to_rank (Mesh_data::NODE);
    face_rank_    = entity_map_->kind_to_rank (Mesh_data::FACE);
    element_rank_ = entity_map_->kind_to_rank (Mesh_data::CELL);


    // Reset all of the mesh-specific data.
    face_id_ = 0;
    Parts (0).swap (element_blocks_);
    Parts (0).swap (side_sets_);
    Parts (0).swap (node_sets_);
    Vector_entity_map ().swap (faces_map_);
    coordinate_field_ = 0;
    set_to_part_.clear ();
    
    // Build the data for the mesh object.

    build_meta_data_ (data, fields);
    build_bulk_data_ (data, cellmap, vertmap, fields);

    Mesh *mesh(NULL);

    // NOTE: Ownership of the pointers meta_data_ and bulk_data_ are
    // transferred to the Mesh instance.  That's where the memory is
    // deleted. IMHO, smart pointers should be used ANYTIME memory is
    // allocated in one code block and deleted somewhere else. I
    // wasted too much of my life looking for memory leaks. -WAP

    mesh = new Mesh (space_dimension, communicator_, entity_map_, 
                     meta_data_, bulk_data_,
                     set_to_part_,
                     *(meta_data_->get_field<Vector_field_type> (std::string ("coordinates"))));
    
    return mesh;
}


void Mesh_factory::build_meta_data_ (const Mesh_data::Data& data, const Mesh_data::Fields& fields)
{
    const int num_element_blocks = data.element_blocks ();
    const int num_side_sets      = data.side_sets ();
    const int num_node_sets      = data.node_sets ();

    const int space_dimension = data.parameters ().dimensions ();

    // Convert element blocks, node and side sets into Parts:

    for (int block = 0; block < num_element_blocks; ++block)
        element_blocks_.push_back (add_element_block_ (data.element_block (block)));

    for (int side_set = 0; side_set < num_side_sets; ++side_set)
        side_sets_.push_back (add_side_set_ (data.side_set (side_set)));

    for (int node_set = 0; node_set < num_node_sets; ++node_set)
        node_sets_.push_back (add_node_set_ (data.node_set (node_set)));

    // Get the universal part. There's only one everything.
    stk::mesh::Part &universal_part (meta_data_->universal_part ());

    // Add a faces part.
    faces_part_ = &meta_data_->declare_part ("Element sides", face_rank_);

    put_coordinate_field_ (universal_part, space_dimension);

    // Declare and Put fields
    for (Mesh_data::Fields::const_iterator field = fields.begin ();
         field != fields.end (); 
         ++field)
    {
        put_field_ (*field, universal_part, space_dimension);
    }

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
void Mesh_factory::build_bulk_data_ (const Mesh_data::Data& data, 
                                     const Epetra_Map& cellmap, 
                                     const Epetra_Map& vertmap,
                                     const Mesh_data::Fields& fields)
{
    const int space_dimension = data.parameters ().dimensions ();

    // declare nodes and elements first, nodes will overlap on processors

    bulk_data_->modification_begin();

    for (int set = 0; set < node_sets_.size (); ++set)
        add_nodes_to_part_ (data.node_set (set), *node_sets_ [set], vertmap);

    for (int block = 0; block < element_blocks_.size (); ++block)
        add_elements_to_part_ (data.element_block (block), *element_blocks_ [block], 
                               cellmap, vertmap);
    
    bulk_data_->modification_end();

    bulk_data_->modification_begin();
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

    std::cerr << "Process " << communicator_.MyPID() 
              << " owns " << nface_local 
              << " of " << nface_global 
              << " faces, index " << global_fidx[communicator_.MyPID()]
              << " to "  << global_fidx[communicator_.MyPID()+1] - 1
              << std::endl;
    ASSERT((global_fidx.back()-1) == nface_global);

    bulk_data_->modification_begin();
    generate_local_faces_(global_fidx[communicator_.MyPID()], false);
    bulk_data_->modification_end();


    // FIXME: Generate cell-to-face relations.  If cell-to-face
    // relations are done here instead of in declare_face_() the same
    // problems happen.  Declaring cell-to-face relations stomps on
    // cell-to-node and/or face-to-node relations.

    // bulk_data_->modification_begin();
    // generate_face_relations_();
    // bulk_data_->modification_end();


    // Put side set faces in the correct parts

    bulk_data_->modification_begin();
    for (int set = 0; set < side_sets_.size (); ++set) {
        add_sides_to_part_ (data.side_set (set), *side_sets_ [set], cellmap);
    }
    bulk_data_->modification_end();

}


void Mesh_factory::put_field_ (const Mesh_data::Field& field_data, 
                               stk::mesh::Part& part, 
                               unsigned int space_dimension)
{
    const unsigned int location = entity_map_->kind_to_rank (field_data.location ());

    if (field_data.type () == Mesh_data::SCALAR)
    {
        Scalar_field_type& field (meta_data_->
                                  declare_field<Scalar_field_type>(field_data.name ()));
        stk::mesh::put_field (field, location, part);
    }

    if (field_data.type () == Mesh_data::VECTOR)
    {
        Vector_field_type& field (meta_data_->
                                  declare_field<Vector_field_type>(field_data.name ()));
        stk::mesh::put_field (field, location, part, space_dimension);
    }

}

void Mesh_factory::put_coordinate_field_ (stk::mesh::Part& part, unsigned int space_dimension)
{
    coordinate_field_ = & meta_data_->declare_field<Vector_field_type>("coordinates");
    stk::mesh::put_field (*coordinate_field_, stk::mesh::Node, part, space_dimension);
}

/** 
 * This routine puts the coordinate fields on @em local nodes.  
 *
 * It's written so that it may be called before or after
 * BulkData::modification_end.
 *
 * The specified @c coordinate_data is for the nodes @em declared on
 * the local processor. When the nodes were declared
 * (::add_nodes_to_part_), shared nodes were declared on multiple
 * processors.  When an entity is declared on multiple processors, the
 * lowest rank processor in that set gets ownership. At least, we can
 * expect that one of the processes in that set got ownership of each
 * shared node.  We can also expect that the list of nodes owned by
 * the local process got all mixed up.  We need to use the specified
 * @c vertmap to identify the current local index given the global
 * index.
 * 
 * @param coordinate_data coordinates of nodes originally declared on this processor
 * @param vertmap local to global node index map for nodes originally declared on this processor
 */
void Mesh_factory::add_coordinates_ (const Mesh_data::Coordinates<double>& coordinate_data,
                                     const Epetra_Map& vertmap)
{

    // Select the local nodes
    stk::mesh::Selector owned(meta_data_->locally_owned_part());
    Entity_vector local_nodes;
    stk::mesh::get_selected_entities (owned,
                                      bulk_data_->buckets (stk::mesh::Node), 
                                      local_nodes);

    // Loop over the local nodes, if the node is owned by this
    // process, set the coordinate
    int node_coordinate_index = 0;
    for (Entity_vector::const_iterator node_it = local_nodes.begin ();
         node_it != local_nodes.end (); ++node_it)
    {
        int global_vidx((*node_it)->identifier());
        ASSERT (vertmap.MyGID(global_vidx));
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


stk::mesh::Part* Mesh_factory::add_element_block_ (const Mesh_data::Element_block& block)
{

    std::ostringstream name;
    if (block.name ().size () > 0)
        name << block.name ();
    else
        name << "element block " << block.id ();

    stk::mesh::Part &new_part (meta_data_->declare_part (name.str (), element_rank_));
    stk::mesh::set_cell_topology 
        (new_part, get_topology_data (block.element_type ()).getTopology ());

    add_set_part_relation_ (block.id (), new_part);

    return &new_part;
}


stk::mesh::Part* Mesh_factory::add_side_set_ (const Mesh_data::Side_set& set)
{
    std::ostringstream name;
    if (set.name ().size () > 0)
        name << set.name ();
    else
        name << "side set " << set.id ();
    stk::mesh::Part &new_part (meta_data_->declare_part (name.str (), face_rank_));

    add_set_part_relation_ (set.id (), new_part);

    return &new_part;
}


stk::mesh::Part* Mesh_factory::add_node_set_ (const Mesh_data::Node_set& set)
{
    std::ostringstream name;
    if (set.name ().size () > 0)
        name << set.name ();
    else
        name << "node set " << set.id ();
    stk::mesh::Part &new_part (meta_data_->declare_part (name.str (), stk::mesh::Node));

    add_set_part_relation_ (set.id (), new_part);


    return &new_part;
}

void Mesh_factory::add_set_part_relation_ (unsigned int set_id, stk::mesh::Part& part)
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
 */
void Mesh_factory::add_elements_to_part_ (const Mesh_data::Element_block& block, stk::mesh::Part &part,
                                          const Epetra_Map& cmap, const Epetra_Map& vmap)
{
    // Add connectivity information via stk::mesh::declare_element
    std::vector<int> storage (block.nodes_per_element ());
    std::vector<int> global_vidx(block.nodes_per_element ());

    for (int local_cidx = 0; local_cidx < block.num_elements (); ++local_cidx)
    {
        block.connectivity (local_cidx, storage.begin ());

        for (unsigned int i = 0; i < block.nodes_per_element (); i++) {
            global_vidx[i] = vmap.GID(storage[i]);
        }

        int global_cidx(cmap.GID(local_cidx));

        try {
            stk::mesh::Entity &element = 
                stk::mesh::declare_element (*bulk_data_, part, global_cidx, &global_vidx[0]);
        } catch (const std::exception& e) {
            std::stringstream msg;
            msg << "cell error: local: " << local_cidx << ": ";
            std::copy(storage.begin(), storage.end(), std::ostream_iterator<int>(msg, ", "));
            msg << "global: " << global_cidx << ": ";
            std::copy(global_vidx.begin(), global_vidx.end(), std::ostream_iterator<int>(msg, ", "));
            std::cerr << msg.str() << std::endl;

            throw e;
        }

    }
}

Entity_vector
Mesh_factory::get_element_side_nodes_(const stk::mesh::Entity& element, 
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
Mesh_factory::get_element_side_face_(const stk::mesh::Entity& element, const unsigned int& s)
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
 * This routine takes care of everything necessary to declare a
 * (local) face, include declaring the entity and making necessary
 * relations to related nodes and cells.  In order for this to work,
 * the elements have to be declared @em and ghosted.  
 * 
 * @param nodes nodes involved in the face
 * @param index global index to assign to the face
 * @param owner the element that is the "owner" of the face
 * @param nbr the neighboring element, if any
 * 
 * @return 
 */
const stk::mesh::Entity&
Mesh_factory::declare_face_(stk::mesh::EntityVector& nodes, const unsigned int& index, 
                            stk::mesh::Entity *owner, const unsigned int& side_index,
                            stk::mesh::Entity *nbr)
{
    stk::mesh::PartVector p;
    p.push_back(faces_part_);

    stk::mesh::Entity& face = 
        bulk_data_->declare_entity(entity_map_->kind_to_rank(Mesh_data::FACE), index, p);

    unsigned int r(0);
    for (stk::mesh::EntityVector::iterator n = nodes.begin(); n != nodes.end(); n++) {
        bulk_data_->declare_relation(face, **n, r++);
    }
    
    // FIXME: declaring cell to face relations here seems to mess up
    // the cell to node relations, not sure why (varying the relation
    // index doesn't seem to matter)

    // bulk_data_->declare_relation(*owner, face, side_index);
    // static const unsigned int ioffset(0);
    
    // bulk_data_->declare_relation(*owner, face, index+ioffset);
    // if (nbr != NULL) {
    //     bulk_data_->declare_relation(*nbr, face, index+ioffset);
    // }

    return face;
}

/** 
 * This routine goes through the locally owned cells and generates
 * cell faces.  The result are faces that are uniquely owned and
 * uniquely defined.  For example, only one face relating two
 * neighboring cells will be declared.
 *
 * Some conventions:
 *
 *  - A face can only be owned by a process that owns the related
 *    cell(s).  
 * 
 *  - If the processes owning cells related to a given face
 *    are different, the face is owned by the lower ranked process.
 *
 * - The order of face-to-cell relations is important.  All faces will
 *    have one relation -- to it's "owner".  The second relation is
 *    for the "neighbor" cell.  This may be used later to determine
 *    face direction.
 *
 * caller is responsible for calling
 * stk::mesh::BulkData::modification_begin() (when @c justcount is @c
 * false ).
 *
 * This would be a lot easier if we could just declare the faces with
 * some bogus global indexes and then renumber them, but I can't
 * figure out how to do that.
 * 
 * @param justcount if true, do not declare any faces, just count how
 * many there would be.
 * 
 * @return number of faces declared/counted
 */
int 
Mesh_factory::generate_local_faces_(const int& faceidx0, const bool& justcount)
{      
    const unsigned int me = communicator_.MyPID();
    int faceidx(faceidx0);
    stk::mesh::Selector owned(meta_data_->locally_owned_part());

    int nlocal = 
        stk::mesh::count_selected_entities(owned, bulk_data_->buckets(entity_map_->kind_to_rank(Mesh_data::CELL)));

    stk::mesh::PartVector parts(meta_data_->get_parts());
    for (stk::mesh::PartVector::iterator p = parts.begin(); 
           p != parts.end(); p++) {
        const CellTopologyData* topo = stk::mesh::get_cell_topology (**p);
        
        if (topo == NULL) continue;

        stk::mesh::Selector s(*(*p));
        s &= owned;

        std::vector< stk::mesh::Entity * > cells;
        stk::mesh::get_selected_entities(s, 
                                         bulk_data_->buckets(entity_map_->kind_to_rank(Mesh_data::CELL)), 
                                         cells);
        
        int localidx(0);
        for (std::vector< stk::mesh::Entity * >::iterator c = cells.begin();
             c != cells.end(); c++, localidx++) {

            int globalidx((*c)->identifier());

            for (unsigned int s = 0; s < topo->side_count; s++) {
                const CellTopologyData *side_topo = (topo->side[s].topology);
                const unsigned * const side_node_map = topo->side[s].node;

                // see if the face already exists (on the local
                // processor); if so, just go on to the next

                stk::mesh::Entity *rface(get_element_side_face_(**c, s));
                if (rface != NULL) continue;
                                         
                // get the cell nodes on this side

                Entity_vector snodes(get_element_side_nodes_(**c, s));

                // get whatever cells the nodes on this side are
                // related to (*c should be one of them, *and* there
                // should be at most one other cell)
            
                Entity_vector rcells;
                get_entities_through_relations(snodes, 
                                               element_rank_, 
                                               rcells);

                if (rcells.size() < 1) {
                    
                    std::string msg = 
                        boost::str(boost::format("error: cell %d (owner = %d): face %d (%d nodes) doesn't even relate to cell") %
                                   (*c)->identifier() % me % s % snodes.size());
                    // std::cerr << msg << std::endl;
                    throw std::runtime_error(msg);
                } else if (rcells.size() > 2) {
                    std::string msg = 
                        boost::str(boost::format("error: cell %d (owner = %d): face %d relates to too many cells (%d)") %
                                   (*c)->identifier() % me % s % rcells.size());
                    throw std::runtime_error(msg);
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
                // rank is less
                
                if (rcell != NULL && rcell->owner_rank() > me) {
                    continue;
                }


                if (!justcount) {
                    const stk::mesh::Entity& face = 
                        declare_face_(snodes, faceidx, *c, s, rcell);
                    // declare_element_side (*bulk_data_, 
                    //                       static_cast<stk::mesh::EntityId>(faceidx),
                    //                       **c, s, faces_part_);
                }
                faceidx++;
                
                std::cerr << "Cell " << globalidx 
                          << ": side " << s 
                          << ": face idx = " << faceidx
                          << std::endl;

            } // side loop

        } // entity (element) loop
 
    } // part loop
    
    return faceidx - faceidx0;
}


void Mesh_factory::add_sides_to_part_ (const Mesh_data::Side_set& side_set, 
                                       stk::mesh::Part &part,
                                       const Epetra_Map& cmap)
{

    // Side set consists of elements (local index) and a local face
    // number. We need to convert these to the unique face indices.

    stk::mesh::PartVector parts_to_add;
    parts_to_add.push_back (&part);

    const int num_sides  = side_set.num_sides ();
    const std::vector<int>& element_list = side_set.element_list ();
    const std::vector<int>& side_list    = side_set.element_list ();
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
            throw std::runtime_error(msg);
        }

        bulk_data_->change_entity_parts(*face, parts_to_add);
    }
}


/** 
 * This routine declares the nodes in the node set and puts them in
 * the specified part.  (NOTE: the nodes should just be declared,
 * don't expect that they were declared elsewhere; entities can be
 * declared many times)
 * 
 * @param node_set 
 * @param part 
 * @param vmap 
 */
void Mesh_factory::add_nodes_to_part_ (const Mesh_data::Node_set& node_set, 
                                       stk::mesh::Part &part,
                                       const Epetra_Map& vmap)
{

    stk::mesh::PartVector parts_to_add;
    parts_to_add.push_back (&part);

    const int num_nodes = node_set.num_nodes ();
    const std::vector<int>& node_list = node_set.node_list ();
    ASSERT (node_list.size () == num_nodes);

    for (std::vector<int>::const_iterator it = node_list.begin ();
         it != node_list.end ();
         ++it)
    {
        int global_vidx(vmap.GID(*it));
        stk::mesh::Entity& node = bulk_data_->declare_entity (stk::mesh::Node, global_vidx, parts_to_add);
    }

}


}


