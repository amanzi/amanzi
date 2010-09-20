#ifndef _MESH_FACTORY_HH_
#define _MESH_FACTORYHH_

#include "Data.hh"
#include "Field_data.hh"
#include "Mesh.hh"

// Trilinos STK_mesh includes.
#include <Shards_BasicTopologies.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/fem/FieldDeclarations.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>

#include <Epetra_MpiComm.h>

namespace Mesh_data
{
class Data;
class Element_block;
class Side_set;
class Node_set;
}

namespace STK_mesh
{

/*!
 * Class which builds an STK_mesh and fields from a Mesh_data::Data
 * and Mesh_data::Fields specification.
 */
class Mesh_factory
{

private:

    const int bucket_size_;
    stk::ParallelMachine parallel_machine_;

    Epetra_MpiComm communicator_;


    typedef std::vector<stk::mesh::Bucket*> Buckets;
    typedef std::vector<stk::mesh::Part*>   Parts;

    typedef std::vector<stk::mesh::EntityId> Entity_id_vector;
    typedef std::map<Entity_id_vector, stk::mesh::Entity*> Vector_entity_map;

    void add_coordinates_ (const Mesh_data::Coordinates<double>& data);
    void build_meta_data_ (const Mesh_data::Data& data, const Mesh_data::Fields& fields);
    void build_bulk_data_ (const Mesh_data::Data& data, const Mesh_data::Fields& fields);
    void receive_bulk_data_ ();

    // Add parts to the meta-data.
    stk::mesh::Part* add_element_block_ (const Mesh_data::Element_block& block);
    stk::mesh::Part* add_side_set_      (const Mesh_data::Side_set& set);
    stk::mesh::Part* add_node_set_      (const Mesh_data::Node_set& set);
    void declare_faces_                 (stk::mesh::Entity& element, stk::mesh::Part &part);

    // Populate parts with elements and fields via the bulk-data
    void add_elements_to_part_ (const Mesh_data::Element_block& block, stk::mesh::Part& part);
    void add_sides_to_part_    (const Mesh_data::Side_set& side_set,   stk::mesh::Part& part);
    void add_nodes_to_part_    (const Mesh_data::Node_set& node_set,   stk::mesh::Part& part);

    void put_field_ (const Mesh_data::Field& field, stk::mesh::Part&, unsigned int space_dimension);

    //! Convert mesh_data field descriptors to STK mesh descriptors.
    stk::mesh::EntityRank map_to_entity_type_ (Mesh_data::Entity_kind location);


    // Temporary information for the mesh currently under construction.

    stk::mesh::BulkData *bulk_data_;
    stk::mesh::MetaData *meta_data_;
    Entity_map          *entity_map_;

    stk::mesh::EntityRank face_rank_;
    stk::mesh::EntityRank element_rank_;
    
    int face_id_;

    Parts element_blocks_;
    Parts side_sets_;
    Parts node_sets_;

    Vector_entity_map faces_map_;
    stk::mesh::Part* faces_part_;

public:

    Mesh_factory (const stk::ParallelMachine& comm, int bucket_size);

    //! Build a mesh from data.
    Mesh* build_mesh (const Mesh_data::Data& data, 
                      const Mesh_data::Fields& fields);

    //! Build a mesh on all the other processors.
    Mesh* build_mesh ();
};

}

#endif
