#ifndef _MESH_MAPS_H_
#define _MESH_MAPS_H_

#include "Mesh.hh"
#include "Entity_map.hh"

#include "dbc.hh"

#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Teuchos_RCPDecl.hpp>

#include <memory>

typedef Teuchos::RCP<STK_mesh::Mesh> Mesh_p;

namespace STK_mesh
{

class Mesh_maps
{

private:

    Mesh_p mesh_;
    const Entity_map& entity_map_;
    Epetra_MpiComm communicator_;

    void update_internals_ ();
    void clear_internals_ ();

    // Maps, Accessors and setters.
    // ----------------------------
    std::auto_ptr<Epetra_Map> maps_ [6];
    const Epetra_Map& map_  (Mesh_data::Entity_kind kind, bool include_ghost) const;
    unsigned int map_index_ (Mesh_data::Entity_kind kind, bool include_ghost) const;
    void assign_map_        (Mesh_data::Entity_kind kind, bool include_ghost, Epetra_Map *map);
    unsigned int kind_to_index_ (Mesh_data::Entity_kind type) const;
    Mesh_data::Entity_kind index_to_kind_ (unsigned int index) const;

    void build_maps_ ();
    void build_tables_ ();

    bool valid_entity_kind_ (Mesh_data::Entity_kind kind) const;

    // Local-id tables of entities
    std::vector<unsigned int> cell_to_face_;
    std::vector<unsigned int> cell_to_node_;
    std::vector<unsigned int> face_to_node_;

    std::map<unsigned int, unsigned int> global_to_local_ [3];

    template <typename F, typename D, typename M>
    void add_global_ids_ (F from, F to, D destination, M& inverse);

public:

    Mesh_maps (Mesh_p mesh);

    void update ();

    // Local id interfaces.
    // --------------------

    template <typename IT>
    void cell_to_faces (unsigned int cell, IT begin, IT end) const;

    template <typename IT>
    void cell_to_nodes (unsigned int cell, IT begin, IT end) const;

    template <typename IT>
    void face_to_nodes (unsigned int face, IT begin, IT end) const;

    template <typename IT>
    void node_to_coordinates (unsigned int node, IT begin, IT end) const;

    template <typename IT>
    void face_to_coordinates (unsigned int face, IT begin, IT end) const;

    template <typename IT>
    void cell_to_ccordinates (unsigned int call, IT begin, IT end) const;

    const Epetra_Map& cell_map (bool include_ghost) const;
    const Epetra_Map& face_map (bool include_ghost) const;
    const Epetra_Map& node_map (bool include_ghost) const;

    unsigned int count_entities (Mesh_data::Entity_kind kind, Element_Category category) const;

    // Entity Sets (cell, side, node)
    // ------------------------------
    unsigned int num_sets (Mesh_data::Entity_kind kind) const;

    template <typename IT>
    void set_ids (Mesh_data::Entity_kind kind, IT begin, IT end) const;

    bool valid_set_id (Mesh_data::Entity_kind kind, unsigned int id) const;

    unsigned int set_size (unsigned int set_id, Mesh_data::Entity_kind kind, 
                           Element_Category category) const;

    template <typename IT>
    void get_set (unsigned int set_id, Mesh_data::Entity_kind kind, Element_Category category,
                  IT begin, IT end) const;


};

// Template & inline members
// ------------------------

const Epetra_Map& Mesh_maps::cell_map (bool include_ghost) const
{
    return map_ (Mesh_data::CELL, include_ghost);
}

const Epetra_Map& Mesh_maps::face_map (bool include_ghost) const
{
    return map_ (Mesh_data::FACE, include_ghost);
}

const Epetra_Map& Mesh_maps::node_map (bool include_ghost) const
{
    return map_ (Mesh_data::NODE, include_ghost);
}


template <typename IT>
void Mesh_maps::cell_to_faces (unsigned int cell, IT destination_begin, IT destination_end) const
{
    ASSERT ((unsigned int) (destination_end - destination_begin) == 6);
    const unsigned int index = 6*cell;
    std::vector<unsigned int>::const_iterator begin = cell_to_face_.begin () + index;
    std::vector<unsigned int>::const_iterator end = begin + 6;
    std::copy (begin, end, destination_begin);
}

template <typename IT>
void Mesh_maps::cell_to_nodes (unsigned int cell, IT destination_begin, IT destination_end) const
{
    ASSERT ((unsigned int) (destination_end - destination_begin) == 8);
    const unsigned int index = 8*cell;
    std::vector<unsigned int>::const_iterator begin = cell_to_node_.begin () + index;
    std::vector<unsigned int>::const_iterator end   = begin + 8;
    std::copy (begin, end, destination_begin);
}

template <typename IT>
void Mesh_maps::face_to_nodes (unsigned int face, IT destination_begin, IT destination_end) const
{
    ASSERT ((unsigned int) (destination_end - destination_begin) == 4);
    const unsigned int index = 4*face;
    std::vector<unsigned int>::const_iterator begin = face_to_node_.begin () + index;
    std::vector<unsigned int>::const_iterator end   = begin + 4;
    std::copy (begin, end, destination_begin);
}


template <typename IT>
void Mesh_maps::node_to_coordinates (unsigned int local_node_id, IT begin, IT end) const
{

    ASSERT ((unsigned int) (end-begin) == 3);

    // Convert local node to global node Id.
    stk::mesh::EntityId global_node_id = map_ (Mesh_data::NODE, true).GID (local_node_id);

    const double * coordinates = mesh_->coordinates (global_node_id);
    std::copy (coordinates, coordinates+3, begin);

}


template <typename F, typename D, typename M>
void Mesh_maps::add_global_ids_ (F from, F to, D destination, M& global_to_local_map)
{
    int local_id = 0;
    for (F it = from; it != to;  ++it)
    {
        const unsigned int global_id = (*it)->identifier ();
        *destination = global_id;
        global_to_local_map.insert (std::make_pair (global_id, local_id));

        destination++;
        local_id++;
    }

}

}



#endif /* _MESH_MAPS_H_ */
