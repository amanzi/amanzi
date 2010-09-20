#include "Mesh_maps.hh"
#include "dbc.hh"

#include <Teuchos_RCP.hpp>
#include <stk_mesh/fem/EntityRanks.hpp>


namespace STK_mesh
{


Mesh_maps::Mesh_maps (Mesh_p mesh) : mesh_ (mesh), 
                                     entity_map_ (mesh->entity_map ()),
                                     communicator_ (mesh_->communicator ())
{
    update();
}

void Mesh_maps::update ()
{
    clear_internals_ ();
    update_internals_ ();
}

void Mesh_maps::clear_internals_ ()
{

    cell_to_face_.resize (0);
    cell_to_node_.resize (0);
    face_to_node_.resize (0);

    for (int i=0; i<3; ++i)
        global_to_local_ [i].erase (global_to_local_ [i].begin (), global_to_local_ [i].end ());

}


void Mesh_maps::update_internals_()
{
    build_maps_ ();
    build_tables_ ();
}


void Mesh_maps::build_maps_ ()
{

    // For each of Elements, Faces, Nodes:
    for (int entity_kind_index = 0; entity_kind_index < 3; ++entity_kind_index)
    {
        
        Mesh_data::Entity_kind kind = index_to_kind_ (entity_kind_index);
        stk::mesh::EntityRank rank = entity_map_.kind_to_rank (kind);

        // Get the collection of "owned" entities 
        Mesh::Entity_vector entities;
        mesh_->get_entities (rank, OWNED, entities);
        const int num_local_entities = entities.size ();
        ASSERT (num_local_entities == mesh_->count_entities (rank, OWNED));

 
        // Get the collection of ghost entities.
        Mesh::Entity_vector ghosts;
        mesh_->get_entities (entity_kind_index, GHOST, ghosts);
        const int num_ghost_entities = ghosts.size ();
        const int num_used_entities = num_local_entities + num_ghost_entities;

        // Put the collections together, ghosts at the back.
        std::copy (ghosts.begin (), ghosts.end (), std::back_inserter (entities));
        ASSERT (entities.size () == num_used_entities);


        // Create a vector of global ids and populate the inverse map.
        std::vector<int> my_entities_global_ids;
        add_global_ids_ (entities.begin (), entities.end (),
                         std::back_inserter (my_entities_global_ids), 
                         global_to_local_ [entity_kind_index]);
        ASSERT (my_entities_global_ids.size () == num_local_entities);
        ASSERT (global_to_local_ [entity_kind_index].size () == num_local_entities);




        // Create the used = local+ghost map.
        Epetra_Map *local_ghost_map (new Epetra_Map (-1,
                                                     num_used_entities,
                                                     &my_entities_global_ids [0],
                                                     0, 
                                                     communicator_));

        assign_map_ (kind, true, local_ghost_map);
        ASSERT (map_ (kind, true).NumMyElements () == num_used_entities);

        // Create the local map.
        // Note that this is not sharing data with the local+ghost map.
        Epetra_Map *local_map (new Epetra_Map (-1,
                                               num_local_entities,
                                               &my_entities_global_ids [0],
                                               0,
                                               communicator_));

        assign_map_ (kind, false, local_map);
        ASSERT (map_ (kind, false).NumMyElements () == num_local_entities);

    }

}

void Mesh_maps::build_tables_ ()
{

    // Cell to faces and nodes.
    ASSERT (cell_to_face_.size () == 0);
    ASSERT (cell_to_node_.size () == 0);

    // Loop over cells using local indices.
    const Epetra_Map &the_cell_map (cell_map (true));
    const unsigned int num_local_cells = count_entities (Mesh_data::CELL, USED);
    for (int local_cell = 0; local_cell < num_local_cells; ++local_cell)
    {
        // Get the global index of the cell from the map.
        const unsigned int global_index = the_cell_map.GID (local_cell);

        Mesh::Entity_Ids faces;
        mesh_->element_to_faces (global_index, faces);
        ASSERT (faces.size () == 6);

        Mesh::Entity_Ids nodes;
        mesh_->element_to_nodes (global_index, nodes);
        ASSERT (nodes.size () == 8);

        // Loop over faces
        for (Mesh::Entity_Ids::const_iterator face = faces.begin ();
             face != faces.end (); ++face)
        {
            const unsigned int face_index = global_to_local_ [0] [*face];
            cell_to_face_.push_back (face_index);
        }

        // Loop over nodes
        for (Mesh::Entity_Ids::const_iterator node = nodes.begin ();
             node != nodes.end (); ++node)
        {
            const unsigned int node_index = global_to_local_ [1] [*node];
            cell_to_node_.push_back (node_index);
        }

    }
    ASSERT (cell_to_face_.size () == 6 * count_entities (Mesh_data::CELL, USED));
    ASSERT (cell_to_node_.size () == 8 * count_entities (Mesh_data::CELL, USED));

    // Faces to nodes
    ASSERT (face_to_node_.size () == 0);
    const Epetra_Map &the_face_map (face_map (true));
    const unsigned int num_local_faces = count_entities (Mesh_data::FACE, USED);
    for (int local_face = 0; local_face < num_local_faces; ++local_face)
    {
        const unsigned int global_index = the_face_map.GID (local_face);
        
        Mesh::Entity_Ids nodes;
        mesh_->face_to_nodes (global_index, nodes);
        ASSERT (nodes.size () == 4);
        
        // Loop over nodes
        for (Mesh::Entity_Ids::const_iterator node = nodes.begin ();
             node != nodes.end (); ++node)
        {
            const unsigned int node_index = global_to_local_ [1] [*node];
            face_to_node_.push_back (node_index);
        }
    }
    ASSERT (face_to_node_.size () == 4 * count_entities (Mesh_data::FACE, USED));

}

const Epetra_Map& Mesh_maps::map_ (Mesh_data::Entity_kind kind, bool include_ghost) const
{
    return *(maps_ [map_index_ (kind, include_ghost)].get ());
}

/* This should accept the entity_kinds for which we have internal data
 */
bool Mesh_maps::valid_entity_kind_ (const Mesh_data::Entity_kind kind) const
{
    return (kind == Mesh_data::NODE) || (kind == Mesh_data::FACE) || (kind == Mesh_data::CELL);
}


unsigned int Mesh_maps::kind_to_index_ (const Mesh_data::Entity_kind kind) const
{
    ASSERT (valid_entity_kind_ (kind));

    if (kind == Mesh_data::NODE) return 0;
    if (kind == Mesh_data::FACE) return 1;
    if (kind == Mesh_data::CELL) return 2;
}

Mesh_data::Entity_kind Mesh_maps::index_to_kind_ (const unsigned int index) const
{
    ASSERT (index >= 0 && index < 3);

    if (index == 0) return Mesh_data::NODE;
    if (index == 1) return Mesh_data::FACE;
    if (index == 2) return Mesh_data::CELL;

}


unsigned int Mesh_maps::map_index_ (Mesh_data::Entity_kind kind, bool include_ghost) const
{
    unsigned int index = kind_to_index_ (kind);
    index = 2*index + (unsigned int) (include_ghost);
    ASSERT (index < 6);
}

void Mesh_maps::assign_map_ (Mesh_data::Entity_kind kind, bool include_ghost, Epetra_Map *map)
{
    const unsigned int map_index = map_index_ (kind, include_ghost);
    maps_ [map_index] = std::auto_ptr<Epetra_Map>(map);
}

unsigned int Mesh_maps::count_entities (Mesh_data::Entity_kind kind, Element_Category category) const
{
    const int rank = entity_map_.kind_to_rank (kind);
    return mesh_->count_entities (rank, category);
}

}
