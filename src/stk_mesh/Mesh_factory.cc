#include "Mesh_factory.hh"
#include "Mesh.hh"

#include "dbc.hh"

#include <iostream>
#include <algorithm>
#include <sstream>
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

    bulk_data_->modification_begin ();

    for (int set = 0; set < node_sets_.size (); ++set)
        add_nodes_to_part_ (data.node_set (set), *node_sets_ [set], vertmap);

    for (int block = 0; block < element_blocks_.size (); ++block)
        add_elements_to_part_ (data.element_block (block), *element_blocks_ [block], 
                               cellmap, vertmap);
    
    // for (int set = 0; set < side_sets_.size (); ++set)
    //     add_sides_to_part_ (data.side_set (set), *side_sets_ [set], cellmap);
        
    bulk_data_->modification_end ();

    // add_coordinates_ (data.coordinates (), vertmap);

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
        double * coordinate_field_data = stk::mesh::field_data (*coordinate_field_, **node_it);
        coordinate_data (local_vidx, coordinate_field_data);
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

void Mesh_factory::declare_faces_ (stk::mesh::Entity& element, stk::mesh::Part &part)
{
    
    const CellTopologyData* element_topology = get_cell_topology (part);
    const int num_faces = element_topology->side_count;

    for (int local_face = 0; local_face < num_faces; ++local_face)
    {

        const CellTopologyData *face_topology     = element_topology->side [local_face].topology;
        unsigned const * const  side_node_map     = element_topology->side [local_face].node;
        const int               face_vertex_count = face_topology->vertex_count;

        // Get the nodes belonging to the element
        stk::mesh::PairIterRelation node_relations = element.relations (stk::mesh::Node);

        // Pick out the nodes for this face.
        std::vector<stk::mesh::EntityId> face_nodes (face_vertex_count);
        for (int vertex = 0; vertex < face_vertex_count; ++vertex)
        {
            const unsigned int local_node = side_node_map [vertex];
            face_nodes [vertex] = node_relations [local_node].entity ()->identifier ();
        }
        std::sort (face_nodes.begin (), face_nodes.end ());
            
        std::pair<Vector_entity_map::iterator, bool> result = 
            faces_map_.insert (std::make_pair (face_nodes, (stk::mesh::Entity*) (NULL)));
        if (result.second)
        {
            ++face_id_;
            stk::mesh::Entity& face = 
                declare_element_side (*bulk_data_, face_id_, element, local_face, faces_part_);
            result.first->second = &face;
        }
        else 
        {
            stk::mesh::Entity* original_face = result.first->second;
            bulk_data_->declare_relation (element, *original_face, local_face);
            faces_map_.erase (result.first);
        }
    }
}


void Mesh_factory::add_sides_to_part_ (const Mesh_data::Side_set& side_set, 
                                       stk::mesh::Part &part,
                                       const Epetra_Map& cmap)
{

    // Side set consists of elements and a local face number. We need
    // to convert these to the unique face indices.

    stk::mesh::PartVector parts_to_add;
    parts_to_add.push_back (&part);

    const int num_sides  = side_set.num_sides ();
    const std::vector<int>& element_list = side_set.element_list ();
    const std::vector<int>& side_list    = side_set.element_list ();
    ASSERT (element_list.size () == num_sides);
    ASSERT (side_list.size ()    == num_sides);

    for (int index=0; index < num_sides; ++index)
    {
        const int element_number = cmap.GID(element_list [index]);
        const int local_side = side_list [index];

        // Look up the element from the Id.
        stk::mesh::Entity *element = bulk_data_->get_entity (element_rank_, element_number);

        // Look up the face from the local face id.
        stk::mesh::PairIterRelation faces = element->relations (face_rank_);
        for (stk::mesh::PairIterRelation::iterator it = faces.begin (); it != faces.end (); ++it)
        {
            if (it->identifier () == local_side)
                bulk_data_->change_entity_parts(*(it->entity ()), parts_to_add);
        }



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


