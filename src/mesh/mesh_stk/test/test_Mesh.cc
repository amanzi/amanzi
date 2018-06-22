#include <UnitTest++.h>

#include "Setup_tests.hh"

#include "Element_block.hh"
#include "Coordinates.hh"
#include "Element_types.hh"
#include "Field_data.hh"

#include <iostream>



SUITE (Mesh)
{

    TEST_FIXTURE (Mesh_setup, Sizes)
    {

        const double *coordinates;

        if (my_pid == 0)
        {
            CHECK_EQUAL (mesh->rank_id (), my_pid);
            CHECK_EQUAL (mesh->count_entities (mesh->kind_to_rank(Amanzi::AmanziMesh::CELL), Amanzi::AmanziMesh::Parallel_type::OWNED), 4);
            CHECK_EQUAL (mesh->count_entities (mesh->kind_to_rank(Amanzi::AmanziMesh::FACE), Amanzi::AmanziMesh::Parallel_type::OWNED), 21);
            CHECK_EQUAL (mesh->count_entities (mesh->kind_to_rank(Amanzi::AmanziMesh::NODE), Amanzi::AmanziMesh::Parallel_type::OWNED), 20);

            CHECK_EQUAL (mesh->count_entities (mesh->kind_to_rank(Amanzi::AmanziMesh::CELL), Amanzi::AmanziMesh::Parallel_type::ALL), 4);
            CHECK_EQUAL (mesh->count_entities (mesh->kind_to_rank(Amanzi::AmanziMesh::FACE), Amanzi::AmanziMesh::Parallel_type::ALL), 21);
            CHECK_EQUAL (mesh->count_entities (mesh->kind_to_rank(Amanzi::AmanziMesh::NODE), Amanzi::AmanziMesh::Parallel_type::ALL), 20);

        }
        else
        {
            CHECK_EQUAL (mesh->rank_id (), my_pid);
            CHECK_EQUAL (mesh->count_entities (mesh->kind_to_rank(Amanzi::AmanziMesh::NODE), Amanzi::AmanziMesh::Parallel_type::OWNED), 0);
            CHECK_EQUAL (mesh->count_entities (mesh->kind_to_rank(Amanzi::AmanziMesh::FACE), Amanzi::AmanziMesh::Parallel_type::OWNED), 0);
            CHECK_EQUAL (mesh->count_entities (mesh->kind_to_rank(Amanzi::AmanziMesh::CELL), Amanzi::AmanziMesh::Parallel_type::OWNED), 0);

            CHECK_EQUAL (mesh->count_entities (mesh->kind_to_rank(Amanzi::AmanziMesh::NODE), Amanzi::AmanziMesh::Parallel_type::ALL), 0);
            CHECK_EQUAL (mesh->count_entities (mesh->kind_to_rank(Amanzi::AmanziMesh::FACE), Amanzi::AmanziMesh::Parallel_type::ALL), 0);
            CHECK_EQUAL (mesh->count_entities (mesh->kind_to_rank(Amanzi::AmanziMesh::CELL), Amanzi::AmanziMesh::Parallel_type::ALL), 0);

        }

        // CHECK_EQUAL (mesh->count_global_entities (mesh->kind_to_rank(Amanzi::AmanziMesh::CELL)), 4);
        // CHECK_EQUAL (mesh->count_global_entities (mesh->kind_to_rank(Amanzi::AmanziMesh::FACE)), 21);
        // CHECK_EQUAL (mesh->count_global_entities (mesh->kind_to_rank(Amanzi::AmanziMesh::NODE)), 20);
    }

    TEST_FIXTURE (Mesh_setup, Coordinates)
    {

        const double* coordinates;

        coordinates = mesh->coordinates (stk::mesh::EntityId (1));
        set_real_coordinates (0,0,0);
        CHECK_ARRAY_EQUAL (coordinates, real_coordinates, 3);

        coordinates = mesh->coordinates (stk::mesh::EntityId (2));
        set_real_coordinates (0, 1, 0);
        CHECK_ARRAY_EQUAL (coordinates, real_coordinates, 3);

        coordinates = mesh->coordinates (stk::mesh::EntityId (3));
        set_real_coordinates (0, 1, 1);
        CHECK_ARRAY_EQUAL (coordinates, real_coordinates, 3);

    
        coordinates = mesh->coordinates (stk::mesh::EntityId (19));
        set_real_coordinates (4.0, 1.0, 1.0);
        CHECK_ARRAY_EQUAL (coordinates, real_coordinates, 3);

        coordinates = mesh->coordinates (stk::mesh::EntityId (20));
        set_real_coordinates (4.0, 0.0, 1.0);
        CHECK_ARRAY_EQUAL (coordinates, real_coordinates, 3);

    }


}
