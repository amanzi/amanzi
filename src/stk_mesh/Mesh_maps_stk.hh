#ifndef _MESH_MAPS_STK_H_
#define _MESH_MAPS_STK_H_

#include "Mesh.hh"
#include "Entity_map.hh"
#include "Data_structures.hh"

#include "dbc.hh"

#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>

#include <memory>

namespace STK_mesh
{

  class Mesh_maps_stk
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

    void build_maps_ ();
    void build_tables_ ();

    stk::mesh::EntityRank kind_to_rank_ (Mesh_data::Entity_kind kind) const {
        return entity_map_.kind_to_rank (kind);
    }

    bool valid_entity_kind_ (Mesh_data::Entity_kind kind) const;

    // Local-id tables of entities
    std::vector<unsigned int> cell_to_face_;
    std::vector<unsigned int> cell_to_node_;
    std::vector<unsigned int> face_to_node_;


    // Global to local index maps and associated bookkeeping.
    Index_map global_to_local_maps_ [3];
    unsigned int kind_to_index_ (Mesh_data::Entity_kind type) const;
    const Index_map& kind_to_map_ (Mesh_data::Entity_kind kind) const;
    Mesh_data::Entity_kind index_to_kind_ (unsigned int index) const;

    unsigned int global_to_local_ (unsigned int global_id, Mesh_data::Entity_kind kind) const;

    // Builds the global->local maps.
    template <typename F, typename D, typename M>
    void add_global_ids_ (F from, F to, D destination, M& inverse);

public:

    Mesh_maps_stk (Mesh_p mesh);

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
    void cell_to_coordinates (unsigned int cell, IT begin, IT end) const;

    inline const Epetra_Map& cell_map (bool include_ghost) const;
    inline const Epetra_Map& face_map (bool include_ghost) const;
    inline const Epetra_Map& node_map (bool include_ghost) const;

    unsigned int count_entities (Mesh_data::Entity_kind kind, Element_Category category) const;

    // Entity Sets (cell, side, node)
    // ------------------------------

    // Number and sizes
    unsigned int num_sets () const;
    unsigned int num_sets (Mesh_data::Entity_kind kind) const;

    unsigned int get_set_size (unsigned int set_id,
                               Mesh_data::Entity_kind kind,
                               Element_Category category) const;

    unsigned int get_set_size (const char* name,
                               Mesh_data::Entity_kind kind,
                               Element_Category category) const;

    // Id numbers
    template <typename IT>
    void get_set_ids (Mesh_data::Entity_kind kind, IT begin, IT end) const;

    bool valid_set_id (unsigned int id, Mesh_data::Entity_kind kind) const;


    template <typename IT>
    void get_set (unsigned int set_id, Mesh_data::Entity_kind kind, Element_Category category,
                  IT begin, IT end) const;

    template <typename IT>
    void get_set (const char* name, Mesh_data::Entity_kind kind, Element_Category category,
                  IT begin, IT end) const;


};

// -------------------------
// Template & inline members
// ------------------------

// Inlined

const Epetra_Map& Mesh_maps_stk::cell_map (bool include_ghost) const
{
    return map_ (Mesh_data::CELL, include_ghost);
}

const Epetra_Map& Mesh_maps_stk::face_map (bool include_ghost) const
{
    return map_ (Mesh_data::FACE, include_ghost);
}

const Epetra_Map& Mesh_maps_stk::node_map (bool include_ghost) const
{
    return map_ (Mesh_data::NODE, include_ghost);
}

// Connectivity accessors
// ----------------------


template <typename IT>
void Mesh_maps_stk::cell_to_faces (unsigned int cell, IT destination_begin, IT destination_end) const
{
    ASSERT ((unsigned int) (destination_end - destination_begin) == 6);
    const unsigned int index = 6*cell;
    std::vector<unsigned int>::const_iterator begin = cell_to_face_.begin () + index;
    std::vector<unsigned int>::const_iterator end = begin + 6;
    std::copy (begin, end, destination_begin);
}

template <typename IT>
void Mesh_maps_stk::cell_to_nodes (unsigned int cell, IT destination_begin, IT destination_end) const
{
    ASSERT ((unsigned int) (destination_end - destination_begin) == 8);
    const unsigned int index = 8*cell;
    std::vector<unsigned int>::const_iterator begin = cell_to_node_.begin () + index;
    std::vector<unsigned int>::const_iterator end   = begin + 8;
    std::copy (begin, end, destination_begin);
}

template <typename IT>
void Mesh_maps_stk::face_to_nodes (unsigned int face, IT destination_begin, IT destination_end) const
{
    ASSERT ((unsigned int) (destination_end - destination_begin) == 4);
    const unsigned int index = 4*face;
    std::vector<unsigned int>::const_iterator begin = face_to_node_.begin () + index;
    std::vector<unsigned int>::const_iterator end   = begin + 4;
    std::copy (begin, end, destination_begin);
}


// Cooordinate Getters
// -------------------

template <typename IT>
void Mesh_maps_stk::node_to_coordinates (unsigned int local_node_id, IT begin, IT end) const
{
    ASSERT ((unsigned int) (end-begin) == 3);

    // Convert local node to global node Id.
    stk::mesh::EntityId global_node_id = map_ (Mesh_data::NODE, true).GID (local_node_id);

    const double * coordinates = mesh_->coordinates (global_node_id);
    std::copy (coordinates, coordinates+3, begin);
}

template <typename IT>
void Mesh_maps_stk::face_to_coordinates (unsigned int local_face_id, IT begin, IT end) const
{
    ASSERT ((unsigned int) (end-begin) == 12);

    unsigned int node_indices [4];
    face_to_nodes (local_face_id, node_indices, node_indices+4);
    for (int i = 0; i < 4; ++i)
    {
        node_to_coordinates (node_indices [i], begin, begin+4);
        begin+=4;
    }


}

template <typename IT>
void Mesh_maps_stk::cell_to_coordinates (unsigned int local_cell_id, IT begin, IT end) const
{
    ASSERT ((unsigned int) (end-begin) == 24);

    unsigned int node_indices [8];
    cell_to_nodes (local_cell_id, node_indices, node_indices+8);
    for (int i = 0; i < 8; ++i)
    {
        node_to_coordinates (node_indices [i], begin, begin+3);
        begin+=3;
    }

}

// Set getters
// -----------

template <typename IT>
void Mesh_maps_stk::get_set_ids (Mesh_data::Entity_kind kind, IT begin, IT end) const
{
    std::vector<unsigned int> ids;
    mesh_->get_set_ids (kind_to_rank_ (kind), ids);
    ASSERT (ids.size () == num_sets (kind));

    IT last = std::copy (ids.begin (), ids.end (), begin);

    ASSERT (last == end);

}

template <typename IT>
void Mesh_maps_stk::get_set (unsigned int set_id, Mesh_data::Entity_kind kind, Element_Category category,
                         IT begin, IT end) const
{
    Entity_vector entities;
    stk::mesh::Part* set_part = mesh_->get_set (set_id, kind_to_rank_ (kind));
    mesh_->get_entities (*set_part, category, entities);
    

    // Convert to local ids.
    for (Entity_vector::const_iterator it = entities.begin ();
         it != entities.end ();
         ++it)
    {
        *begin = global_to_local_ ( (*it)->identifier (), kind);
        begin++;
    }

    ASSERT (begin == end);
}

template <typename IT>
void Mesh_maps_stk::get_set (const char* name, Mesh_data::Entity_kind kind, Element_Category category,
                         IT begin, IT end) const
{

    Entity_vector entities;
    stk::mesh::Part* set_part = mesh_->get_set (name, kind_to_rank_ (kind));
    mesh_->get_entities (*set_part, category, entities);
    
    // Convert to local ids.
    for (Entity_vector::const_iterator it = entities.begin ();
         it != entities.end ();
         ++it)
    {
        *begin = global_to_local_ ( (*it)->identifier (), kind);
        begin++;
    }

    ASSERT (begin == end);
    

}


// Internal template functions
// ---------------------------


template <typename F, typename D, typename M>
void Mesh_maps_stk::add_global_ids_ (F from, F to, D destination, M& global_to_local_map)
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
