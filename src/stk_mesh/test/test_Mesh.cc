#include <UnitTest++.h>

#include "Setup_tests.hh"

#include "../Mesh.hh"
#include "../Mesh_factory.hh"
#include "../Mesh_maps_stk.hh"
#include "../Element_category.hh"

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

        // CHECK_EQUAL (mesh->count_global_entities (stk::mesh::Element), 4);
        // CHECK_EQUAL (mesh->count_global_entities (stk::mesh::Face), 21);
        // CHECK_EQUAL (mesh->count_global_entities (stk::mesh::Node), 20);

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
