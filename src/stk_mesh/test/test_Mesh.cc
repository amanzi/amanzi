#include <UnitTest++.h>

#include "Test_Mesh.hh"

#include "../Mesh.hh"
#include "../Mesh_factory.hh"
#include "../Mesh_maps.hh"
#include "../Element_category.hh"

#include "Element_block.hh"
#include "Coordinates.hh"
#include "Element_types.hh"
#include "Field_data.hh"

#include <stk_util/parallel/Parallel.hpp>
#include <Teuchos_RCP.hpp>

#include <Epetra_MpiComm.h>

#include <iostream>

typedef Teuchos::RCP<STK_mesh::Mesh> Mesh_p;
typedef Teuchos::RCP<STK_mesh::Mesh_factory> Factory_p;
typedef Teuchos::RCP<Epetra_MpiComm> Comm_p;

struct Mesh_setup : public Test_mesh
{

    Factory_p factory;
    Mesh_p mesh;
    Comm_p communicator;
    int my_pid;
    

    Mesh_setup ()
    {
        Mesh_data::Fields fields;

        int argc = 1;
        char **argv;
        argv = new char*[2];
        argv [0] = (char*)"random_dent.exe";
        argv [1] = NULL;

        stk::ParallelMachine parallel_machine = stk::parallel_machine_init (&argc, &argv);
        communicator = Comm_p (new Epetra_MpiComm(parallel_machine));
        my_pid = communicator->MyPID ();

        const int bucket_size = 20;
        
        factory = Factory_p(new STK_mesh::Mesh_factory (parallel_machine, bucket_size));
        mesh    = Mesh_p (factory->build_mesh (*data, fields));

    }

    ~Mesh_setup ()
    {
        communicator->Barrier ();
        stk::parallel_machine_finalize ();
    }

};



TEST_FIXTURE (Mesh_setup, Build)
{

    if (my_pid == 0)
    {
        CHECK_EQUAL (mesh->rank_id (), my_pid);
        CHECK_EQUAL (mesh->count_entities (stk::mesh::Element, STK_mesh::OWNED), 4);
        CHECK_EQUAL (mesh->count_entities (stk::mesh::Face,    STK_mesh::OWNED), 21);
        CHECK_EQUAL (mesh->count_entities (stk::mesh::Node,    STK_mesh::OWNED), 20);

        CHECK_EQUAL (mesh->count_entities (stk::mesh::Element, STK_mesh::USED), 4);
        CHECK_EQUAL (mesh->count_entities (stk::mesh::Face,    STK_mesh::USED), 21);
        CHECK_EQUAL (mesh->count_entities (stk::mesh::Node,    STK_mesh::USED), 20);

    }
    else
    {
        CHECK_EQUAL (mesh->rank_id (), my_pid);
        CHECK_EQUAL (mesh->count_entities (stk::mesh::Node,    STK_mesh::OWNED), 0);
        CHECK_EQUAL (mesh->count_entities (stk::mesh::Face,    STK_mesh::OWNED), 0);
        CHECK_EQUAL (mesh->count_entities (stk::mesh::Element, STK_mesh::OWNED), 0);

        CHECK_EQUAL (mesh->count_entities (stk::mesh::Node,    STK_mesh::USED), 0);
        CHECK_EQUAL (mesh->count_entities (stk::mesh::Face,    STK_mesh::USED), 0);
        CHECK_EQUAL (mesh->count_entities (stk::mesh::Element, STK_mesh::USED), 0);

    }

    std::vector<unsigned int> faces (6);
    std::vector<unsigned int> nodes (8);
    std::vector<unsigned int> face_nodes (4);

    // These result arrays assume that the global ordering is
    // preserved when assigning local ids and respects the relevant
    // shards ordering.
    int faces_cell_0 [] = {0, 1, 2, 3, 4, 5};
    int nodes_cell_0 [] = {0, 4, 5, 1, 3, 7, 6, 2};

    int nodes_face_0 [] = {0, 4, 7, 3};
        
    STK_mesh::Mesh_maps mesh_map(mesh);
    
    if (my_pid == 0)
    {
    
        mesh_map.cell_to_faces (0, faces.begin (), faces.end ());
        CHECK_ARRAY_EQUAL (faces.begin (), faces_cell_0, 6);

        mesh_map.cell_to_nodes (0, nodes.begin (), nodes.end ());
        CHECK_ARRAY_EQUAL (nodes.begin (), nodes_cell_0, 8);

        mesh_map.face_to_nodes (0, face_nodes.begin (), face_nodes.end ());
        CHECK_ARRAY_EQUAL (face_nodes.begin (), nodes_face_0, 4);

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
