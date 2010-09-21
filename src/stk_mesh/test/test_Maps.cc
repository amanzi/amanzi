#include <UnitTest++.h>

#include "Setup_tests.hh"

SUITE (MAPS)
{

    TEST_FIXTURE (Map_setup, Maps)
    {

        if (my_pid == 0)
        {
    
            const Epetra_Map &local_node_map    (mesh_map.node_map (false));
            const Epetra_Map &local_face_map    (mesh_map.face_map (false));
            const Epetra_Map &local_element_map (mesh_map.cell_map (false));

            CHECK_EQUAL (local_node_map.NumMyElements    (), 20);
            CHECK_EQUAL (local_face_map.NumMyElements    (), 21);
            CHECK_EQUAL (local_element_map.NumMyElements (), 4);
        
            const Epetra_Map &complete_node_map    (mesh_map.node_map (true));
            const Epetra_Map &complete_face_map    (mesh_map.face_map (true));
            const Epetra_Map &complete_element_map (mesh_map.cell_map (true));

            CHECK_EQUAL (complete_node_map.NumMyElements    (), 20);
            CHECK_EQUAL (complete_face_map.NumMyElements    (), 21);
            CHECK_EQUAL (complete_element_map.NumMyElements (), 4);

            CHECK_EQUAL (mesh_map.count_entities (Mesh_data::NODE, STK_mesh::OWNED), 20);
            CHECK_EQUAL (mesh_map.count_entities (Mesh_data::FACE, STK_mesh::OWNED), 21);
            CHECK_EQUAL (mesh_map.count_entities (Mesh_data::CELL, STK_mesh::OWNED), 4);

        }
    }

    TEST_FIXTURE (Map_setup, Coordinates)
    {
        if (my_pid == 0)
        {

            double coord_array [3];

            mesh_map.node_to_coordinates (1, coord_array, coord_array+3);
            set_real_coordinates (0.0, 0.0, 0.0);
            CHECK_ARRAY_EQUAL (coord_array, real_coordinates, 3);

            mesh_map.node_to_coordinates (19, coord_array, coord_array+3);
            set_real_coordinates (4.0, 1.0, 1.0);
            CHECK_ARRAY_EQUAL (coord_array, real_coordinates, 3);

        }
    }

    TEST_FIXTURE (Map_setup, Entity_Relations)
    {

        // These result arrays assume that the global ordering is
        // preserved when assigning local ids and respects the relevant
        // shards ordering.
        int faces_cell_0 [] = {0, 1, 2, 3, 4, 5};
        int nodes_cell_0 [] = {0, 4, 5, 1, 3, 7, 6, 2};

        int nodes_face_0 [] = {0, 4, 7, 3};
        

        std::vector<unsigned int> faces (6);
        std::vector<unsigned int> nodes (8);
        std::vector<unsigned int> face_nodes (4);


        if (my_pid == 0)
        {
            mesh_map.cell_to_faces (0, faces.begin (), faces.end ());
            CHECK_ARRAY_EQUAL (faces.begin (), faces_cell_0, 6);

            mesh_map.cell_to_nodes (0, nodes.begin (), nodes.end ());
            CHECK_ARRAY_EQUAL (nodes.begin (), nodes_cell_0, 8);

            mesh_map.face_to_nodes (0, face_nodes.begin (), face_nodes.end ());
            CHECK_ARRAY_EQUAL (face_nodes.begin (), nodes_face_0, 4);


        }

    }

}
