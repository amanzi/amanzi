#include "Mesh_factory.hh"
#include "Mesh.hh"
#include "Mesh_maps_stk.hh"
#include "Exodus_readers.hh"
#include "Data.hh"

#include <stk_util/parallel/Parallel.hpp>



Mesh_data::Data* read_file (const char* filename)
{

    Mesh_data::Data *data = ExodusII::read_exodus_file (filename);

    return data;
}

void dump_mesh (const STK_mesh::Mesh& mesh)
{

    std::cout << "Mesh dump:" << std::endl;
    std::cout << "----------" << std::endl;
    std::cout << "Number of nodes:    " << mesh.count_entities (stk::mesh::Node,    OWNED) << std::endl;
    std::cout << "Number of faces:    " << mesh.count_entities (stk::mesh::Face,    OWNED) << std::endl;
    std::cout << "Number of elements: " << mesh.count_entities (stk::mesh::Element, OWNED) << std::endl;
    std::cout << std::endl;

    std::cout << "Number of sets: " << std::endl;
    std::cout << "Total: "    << mesh.num_sets ()                   << std::endl;
    std::cout << "Nodes: "    << mesh.num_sets (stk::mesh::Node)    << std::endl;
    std::cout << "Faces: "    << mesh.num_sets (stk::mesh::Face)    << std::endl;
    std::cout << "Elements: " << mesh.num_sets (stk::mesh::Element) << std::endl;
    std::cout << std::endl;

    std::cout << "Set info:" << std::endl;
    for (STK_mesh::Id_map::const_iterator it = mesh.sets_begin ();
         it != mesh.sets_end ();
         ++it)
    {

        const STK_mesh::Rank_and_id set_id = it->first;
        const stk::mesh::Part*      part   = it->second;
        std::cout << "  Kind: "        << set_id.first
                  << "    Id: "        << set_id.second
                  << "    Part name: " << part->name ()
                  << "    Size: "      << mesh.count_entities (*part, OWNED)
                  << std::endl;

    }

    std::cout << std::endl << std::endl;


}

void dump_maps (const STK_mesh::Mesh_maps_stk& maps)
{

    std::cout << "Maps dump:" << std::endl;
    std::cout << "----------" << std::endl;

    std::cout << "Number of nodes: " << maps.count_entities (Mesh_data::NODE, OWNED) << std::endl;
    std::cout << "Number of faces: " << maps.count_entities (Mesh_data::FACE, OWNED) << std::endl;
    std::cout << "Number of cells: " << maps.count_entities (Mesh_data::CELL, OWNED) << std::endl;
    std::cout << std::endl;

    std::cout << "Number of sets: " << std::endl;
    std::cout << "Total: " << maps.num_sets ()                << std::endl;
    std::cout << "Nodes: " << maps.num_sets (Mesh_data::NODE) << std::endl;
    std::cout << "Faces: " << maps.num_sets (Mesh_data::FACE) << std::endl;
    std::cout << "Cells: " << maps.num_sets (Mesh_data::CELL) << std::endl;
    std::cout << std::endl;

    std::cout << "Element blocks: " << std::endl;

    const unsigned int num_element_blocks = maps.num_sets (Mesh_data::CELL);
    std::cout << "  Count: " << num_element_blocks << std::endl;
    std::vector<unsigned int> block_ids (num_element_blocks);

    maps.get_set_ids (Mesh_data::CELL, block_ids.begin (), block_ids.end ());

    for (std::vector<unsigned int>::const_iterator it = block_ids.begin ();
         it != block_ids.end ();
         ++it)
    {
        const unsigned int block_size = maps.get_set_size (*it, Mesh_data::CELL, OWNED);
        std::vector<unsigned int> local_cell_numbers (block_size);
        maps.get_set (*it, Mesh_data::CELL, OWNED, local_cell_numbers.begin (), local_cell_numbers.end ());

        std::ostringstream word_up;
        word_up << "Element block: " << (*it) << "\n"
                << "  Size: " << block_size << "\n"
                << "  Cell ids: ";
        std::copy (local_cell_numbers.begin (), local_cell_numbers.end (),
                   std::ostream_iterator<unsigned int>(word_up, ", "));

        std::cout << word_up.str () << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Side sets:" << std::endl;
    const unsigned int num_side_sets  = maps.num_sets (Mesh_data::FACE);
    std::cout << "  Count: " << num_side_sets << std::endl;
    std::vector<unsigned int> side_set_ids (num_side_sets);
    
    maps.get_set_ids (Mesh_data::FACE, side_set_ids.begin (), side_set_ids.end ());
    
    for (std::vector<unsigned int>::const_iterator it = side_set_ids.begin ();
         it != side_set_ids.end ();
         ++it)
    {
        const unsigned int side_set_size = maps.get_set_size (*it, Mesh_data::FACE, OWNED);
        std::vector<unsigned int> local_face_numbers (side_set_size);
        maps.get_set (*it, Mesh_data::FACE, OWNED, local_face_numbers.begin (), local_face_numbers.end ());

        std::ostringstream word_up;
        word_up << "Side set: " << (*it) << "\n"
                << "  Size: "   << side_set_size << "\n";
        std::cout << word_up.str () << std::endl;
        
    }


}


int main(int argc, char *argv[])
{

    stk::ParallelMachine parallel_machine = stk::parallel_machine_init (&argc, &argv);
    const int bucket_size = 100;

    STK_mesh::Mesh_factory factory (parallel_machine, bucket_size);
    Mesh_data::Fields fields;

    for (int i = 1; i < argc; ++i)
    {
        Mesh_data::Data *data = read_file (argv [i]);
        STK_mesh::Mesh_p mesh = STK_mesh::Mesh_p (factory.build_mesh (*data, fields));
        STK_mesh::Mesh_maps_stk mesh_map (mesh);

        dump_mesh (*mesh);
        dump_maps (mesh_map);
    }

    stk::parallel_machine_finalize ();

    return 0;
}
