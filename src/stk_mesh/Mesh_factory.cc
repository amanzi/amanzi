#include "Mesh_factory.hh"
#include "Mesh.hh"

#include "dbc.hh"

#include <sstream>

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


Mesh* Mesh_factory::build_mesh (const Mesh_data::Data& data, 
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
    Parts (0).swap     (element_blocks_);
    Parts (0).swap     (side_sets_);
    Parts (0).swap     (node_sets_);
    Vector_entity_map ().swap (faces_map_);
    coordinate_field_ = 0;

    build_meta_data_ (data, fields);
    build_bulk_data_ (data, fields);

    Mesh *mesh = new Mesh (space_dimension, communicator_, entity_map_, meta_data_, bulk_data_,
                           *(meta_data_->get_field<Mesh::Vector_field_type> (std::string ("coordinates"))));

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

void Mesh_factory::build_bulk_data_ (const Mesh_data::Data& data, const Mesh_data::Fields& fields)
{
    // Now add contents of the blocks and sets to their respctive parts:

    const int space_dimension = data.parameters ().dimensions ();

    bulk_data_->modification_begin ();

    // Explicitly assuming that all of the mesh data is on node 0:
    if (communicator_.MyPID () == 0)
    {

        for (int block = 0; block < element_blocks_.size (); ++block)
            add_elements_to_part_ (data.element_block (block), *element_blocks_ [block]);
        
        for (int set = 0; set < side_sets_.size (); ++set)
            add_sides_to_part_ (data.side_set (set), *side_sets_ [set]);
        
        for (int set = 0; set < node_sets_.size (); ++set)
            add_nodes_to_part_ (data.node_set (set), *node_sets_ [set]);
    }

    bulk_data_->modification_end ();

    add_coordinates_ (data.coordinates ());

    // Loop over fields. Add them.
        

}


void Mesh_factory::put_field_ (const Mesh_data::Field& field_data, 
                               stk::mesh::Part& part, 
                               unsigned int space_dimension)
{
    const unsigned int location = entity_map_->kind_to_rank (field_data.location ());

    if (field_data.type () == Mesh_data::SCALAR)
    {
        Mesh::Scalar_field_type& field (meta_data_->
                                        declare_field<Mesh::Scalar_field_type>(field_data.name ()));
        stk::mesh::put_field (field, location, part);
    }

    if (field_data.type () == Mesh_data::VECTOR)
    {
        Mesh::Vector_field_type& field (meta_data_->
                                        declare_field<Mesh::Vector_field_type>(field_data.name ()));
        stk::mesh::put_field (field, location, part, space_dimension);
    }

}

void Mesh_factory::put_coordinate_field_ (stk::mesh::Part& part, unsigned int space_dimension)
{
    coordinate_field_ = & meta_data_->declare_field<Mesh::Vector_field_type>("coordinates");
    stk::mesh::put_field (*coordinate_field_, stk::mesh::Node, part, space_dimension);
}


void Mesh_factory::add_coordinates_ (const Mesh_data::Coordinates<double>& coordinate_data)
{

    const int nodes = coordinate_data.nodes ();

    // Select all of the nodes. 
    stk::mesh::Selector universal_selector (meta_data_->universal_part ());
    Entity_vector all_nodes;
    stk::mesh::get_selected_entities (universal_selector, bulk_data_->buckets (stk::mesh::Node), all_nodes);
    ASSERT (all_nodes.size () == nodes);

    // Loop over all nodes
    int node_coordinate_index = 0;
    for (Entity_vector::const_iterator node_it = all_nodes.begin ();
         node_it != all_nodes.end ();
         ++node_it)
    {
        double * coordinate_field_data = stk::mesh::field_data (*coordinate_field_, **node_it);

        coordinate_data (node_coordinate_index, coordinate_field_data);
        node_coordinate_index++;
    }

}


stk::mesh::Part* Mesh_factory::add_element_block_ (const Mesh_data::Element_block& block)
{
    stk::mesh::Part &new_part (meta_data_->declare_part (block.name (), element_rank_));
    stk::mesh::set_cell_topology 
        (new_part, get_topology_data (block.element_type ()).getTopology ());
    return &new_part;
}


stk::mesh::Part* Mesh_factory::add_side_set_ (const Mesh_data::Side_set& set)
{
    std::ostringstream name ("Side set: ");
    name << set.id ();
    stk::mesh::Part &new_part (meta_data_->declare_part (name.str (), face_rank_));
    return &new_part;
}


stk::mesh::Part* Mesh_factory::add_node_set_ (const Mesh_data::Node_set& set)
{
    std::ostringstream name ("Node set: ");
    name << set.id ();
    stk::mesh::Part &new_part (meta_data_->declare_part (name.str (), stk::mesh::Node));
    return &new_part;
}


void Mesh_factory::add_elements_to_part_ (const Mesh_data::Element_block& block, stk::mesh::Part &part)
{
    // Add connectivity information via stk::mesh::declare_element
    std::vector<int> storage (block.nodes_per_element ());
    for (int element_ind = 0; element_ind < block.num_elements (); ++element_ind)
    {
        block.connectivity (element_ind, storage.begin ());
        stk::mesh::Entity &element = 
            stk::mesh::declare_element (*bulk_data_, part, (element_ind+1), &storage [0]);

        // XXX We don't need this for general element blocks.
        declare_faces_ (element, part);
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


void Mesh_factory::add_sides_to_part_ (const Mesh_data::Side_set& side_set, stk::mesh::Part &part)
{

}


void Mesh_factory::add_nodes_to_part_ (const Mesh_data::Node_set& node_set, stk::mesh::Part &part)
{

}


}
