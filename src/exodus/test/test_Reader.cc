#include <UnitTest++.h>

#include <stdexcept>
#include <memory>
#include <errno.h>

#include "../Exodus_file.hh"
#include "../Exodus_readers.hh"

#include "Parameters.hh"
#include "Side_set.hh"

struct Exodus_file_holder
{
    ExodusII::Exodus_file file;
    std::auto_ptr<Mesh_data::Data> data;

    Exodus_file_holder (const char* filename) :
        file (filename),
        data (ExodusII::read_exodus_file (filename))
    {  }
};

struct Big_File : Exodus_file_holder
{
    Big_File () : Exodus_file_holder ("exodus/test_files/htc_rad_test-random.exo") { }
};

struct quad_4x4 : Exodus_file_holder
{
    quad_4x4 () : Exodus_file_holder ("exodus/test_files/quad_4x4_ss.exo") { }
};


SUITE (Big_File)
{

    TEST_FIXTURE (Big_File, Parameters)
    {
        const Mesh_data::Parameters &params (data->parameters ());

        CHECK_EQUAL (params.element_block_ids_.size (), 3);
        CHECK_EQUAL (params.node_set_ids_.size (), 0);
        CHECK_EQUAL (params.side_set_ids_.size (), 4);

        CHECK_EQUAL (params.dimensions_, 3);
        CHECK_EQUAL (params.num_nodes_, 6615);
        CHECK_EQUAL (params.num_elements_, 5600);
        CHECK_EQUAL (params.num_element_blocks_, 3);
        CHECK_EQUAL (params.num_node_sets_, 0);
        CHECK_EQUAL (params.num_side_sets_, 4);
    };

    TEST_FIXTURE (Big_File, Side_Sets)
    {

        {
            const Mesh_data::Side_set &side (data->side_set (0));

            // CHECK       (side.has_node_factors ());
            CHECK_EQUAL (side.num_sides (), 400);
            // CHECK_EQUAL (side.num_nodes (), 1600);
        }

        {
            const Mesh_data::Side_set &side (data->side_set (1));

            // CHECK       (side.has_node_factors ());
            CHECK_EQUAL (side.num_sides (), 400);
            // CHECK_EQUAL (side.num_nodes (), 1600);
        }

        {
            const Mesh_data::Side_set &side (data->side_set (2));

            // CHECK       (side.has_node_factors ());
            CHECK_EQUAL (side.num_sides (), 1120);
            // CHECK_EQUAL (side.num_nodes (), 4480);
        }

        {
            const Mesh_data::Side_set &side (data->side_set (3));

            // CHECK       (side.has_node_factors ());
            CHECK_EQUAL (side.num_sides (), 800);
            // CHECK_EQUAL (side.num_nodes (), 3200);
        }

    }

}

SUITE (quad_4x4)
{

    TEST_FIXTURE (quad_4x4, Parameters)
    {
        const Mesh_data::Parameters &params (data->parameters ());

        CHECK_EQUAL (params.element_block_ids_.size (), 3);
        CHECK_EQUAL (params.node_set_ids_.size (), 0);
        CHECK_EQUAL (params.side_set_ids_.size (), 4);

        CHECK_EQUAL (params.dimensions_, 3);
        CHECK_EQUAL (params.num_nodes_, 6615);
        CHECK_EQUAL (params.num_elements_, 5600);
        CHECK_EQUAL (params.num_element_blocks_, 3);
        CHECK_EQUAL (params.num_node_sets_, 0);
        CHECK_EQUAL (params.num_side_sets_, 4);
    };

    TEST_FIXTURE (quad_4x4, Side_Sets)
    {

        {
            const Mesh_data::Side_set &side (data->side_set (0));

            // CHECK       (side.has_node_factors ());
            CHECK_EQUAL (side.num_sides (), 400);
            // CHECK_EQUAL (side.num_nodes (), 1600);
        }

        {
            const Mesh_data::Side_set &side (data->side_set (1));

            // CHECK       (side.has_node_factors ());
            CHECK_EQUAL (side.num_sides (), 400);
            // CHECK_EQUAL (side.num_nodes (), 1600);
        }

        {
            const Mesh_data::Side_set &side (data->side_set (2));

            // CHECK       (side.has_node_factors ());
            CHECK_EQUAL (side.num_sides (), 1120);
            // CHECK_EQUAL (side.num_nodes (), 4480);
        }

        {
            const Mesh_data::Side_set &side (data->side_set (3));

            // CHECK       (side.has_node_factors ());
            CHECK_EQUAL (side.num_sides (), 800);
            // CHECK_EQUAL (side.num_nodes (), 3200);
        }

    }



}



