#include <UnitTest++.h>

#include <stdexcept>
#include <memory>
#include <errno.h>

#include "../Exodus_file.hh"
#include "../Exodus_readers.hh"

#include "Parameters.hh"
#include "Side_set.hh"

struct Big_File
{
    
    ExodusII::Exodus_file file;
    std::auto_ptr<Mesh_data::Parameters> params;
    
    Big_File () : file ("exodus/test_files/htc_rad_test-random.exo"), params (read_parameters (file))
    { 
        
    }

};

SUITE (Params)
{

    TEST_FIXTURE (Big_File, Sizes)
    {
        CHECK_EQUAL (params->dimensions_, 3);
        CHECK_EQUAL (params->num_nodes_, 6615);
        CHECK_EQUAL (params->num_elements_, 5600);
        CHECK_EQUAL (params->num_element_blocks_, 3);
        CHECK_EQUAL (params->num_node_sets_, 0);
        CHECK_EQUAL (params->num_side_sets_, 4);
    };

    TEST_FIXTURE (Big_File, Data)
    {
        
        CHECK_EQUAL (params->element_block_ids_.size (), 3);
        CHECK_EQUAL (params->node_set_ids_.size (), 0);
        CHECK_EQUAL (params->side_set_ids_.size (), 4);

    }

}



SUITE (Side_set)
{

    TEST_FIXTURE (Big_File, Side_1)
    {
        std::auto_ptr<Mesh_data::Side_set> side ( read_side_set (file, 1));
        
        CHECK       (side->has_node_factors ());
        CHECK_EQUAL (side->num_sides (), 400);
        CHECK_EQUAL (side->num_nodes (), 1600);

    }

    TEST_FIXTURE (Big_File, Side_2)
    {
        std::auto_ptr<Mesh_data::Side_set> side ( read_side_set (file, 2));
        
        CHECK       (side->has_node_factors ());
        CHECK_EQUAL (side->num_sides (), 400);
        CHECK_EQUAL (side->num_nodes (), 1600);

    }

    TEST_FIXTURE (Big_File, Side_3)
    {
        std::auto_ptr<Mesh_data::Side_set> side ( read_side_set (file, 3));
        
        CHECK       (side->has_node_factors ());
        CHECK_EQUAL (side->num_sides (), 1120);
        CHECK_EQUAL (side->num_nodes (), 4480);

    }

    TEST_FIXTURE (Big_File, Side_4)
    {
        std::auto_ptr<Mesh_data::Side_set> side ( read_side_set (file, 4));
        
        CHECK       (side->has_node_factors ());
        CHECK_EQUAL (side->num_sides (), 800);
        CHECK_EQUAL (side->num_nodes (), 3200);

    }

}



SUITE (Element_block)
{

    
}

SUITE (Entire_mesh)
{

    TEST(Data)
    {
        
        Mesh_data::Data* data = ExodusII::read_exodus_file ("htc_rad_test-random.exo");

    }



}
