#include "Mesh_factory.hh"
#include "Mesh.hh"
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

    std::cout << "Number of elements: " << mesh.count_entities (stk::mesh::Element, STK_mesh::OWNED) << std::endl;
    std::cout << "Number of faces:    " << mesh.count_entities (stk::mesh::Face,    STK_mesh::OWNED) << std::endl;
    std::cout << "Number of nodes:    " << mesh.count_entities (stk::mesh::Node,    STK_mesh::OWNED) << std::endl;
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

        const stk::mesh::Part* part = it->second;
        const STK_mesh::Rank_and_id set_id = it->first;
        std::cout << "Kind: " <<  set_id.first << "  Id: " << set_id.second << "  Part name: " << part->name ()
                  << std::endl;

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
        STK_mesh::Mesh *mesh = factory.build_mesh (*data, fields);
        
        dump_mesh (*mesh);
        delete mesh;

    }

    stk::parallel_machine_finalize ();


    return 0;
}
