#ifndef _EXAMPLE_MESH_H_
#define _EXAMPLE_MESH_H_


#include "Element_block.hh"
#include "Coordinates.hh"
#include "Element_types.hh"
#include "Data.hh"

/* Global node ordering
 *
 *      3       7      11      15      19         
 *      +-------+-------+-------+-------+        
 *     /       /       /       /       /|       
 *   4/      8/     12/     16/     20/ |       
 *   +-------+-------+-------+-------+  |       
 *   |       |       |       |       |  +18        Z  Y
 *   |  e1   |  e2   |  e3   |  e4   | /           | /
 *   |       |       |       |       |/            |/
 *   +-------+-------+-------+-------+             *--X
 *   1       5      9       13      17          
 *
 *  Local node numbering
 *      8       7
 *      +-------+
 *     /       /| 
 *   5/      6/ | 
 *   +-------+  |
 *   |       |  +3 
 *   |  e1   | /
 *   |       |/ 
 *   +-------+
 *   1       2  
 */



const int local_global_map [] = {1, 5, 6, 2, 4, 8, 7, 3};
static void make_node_ids (int element_number, int node_ids [])
{
    const int base = 4*element_number;
    for (int node = 0; node < 8; ++node)
        node_ids [node] = base+local_global_map [node];
}

struct Test_mesh
{

    Mesh_data::Data *data;


    Test_mesh ()
    {
        const int element_block_id = 1;

        const int dimensions = 3;
        const int num_elements = 4;
        const int num_nodes = 20;
        const int num_blocks = 1;
        const int num_side_sets = 0;
        const int num_node_sets = 0;

        // Parameters
        Mesh_data::Parameters *parameters = 
            new Mesh_data::Parameters("Test mesh", 
                                      dimensions, 
                                      num_nodes, 
                                      num_elements, 
                                      num_blocks, num_side_sets, num_node_sets, 
                                      std::vector<int>(1,element_block_id),
                                      std::vector<int>(0),
                                      std::vector<int>(0));


        // Element block
        std::vector<int> connectivity (32);
        for (int element=0; element<4; ++element)
        {
            make_node_ids (element, &connectivity [element*8]);
        }
        std::vector<double> attributes(0);
        std::string name ("A");
        Mesh_data::Element_block* block = 
            Mesh_data::Element_block::build_from(element_block_id, 
                                                 name,
                                                 num_elements, 
                                                 Mesh_data::HEX, 
                                                 connectivity, attributes);

        std::vector<Mesh_data::Element_block*> element_blocks(1, block);



        // Coordinates
        std::vector<std::vector<double> > coord_data (3, std::vector<double>(20));
        for (int node = 0; node < 20; ++node)
        {
            const int x_plane  = node / 4;
            const int in_plane = node % 4;
            const int y_plane = int ((in_plane == 1) || (in_plane == 2));
            const int z_plane = int ((in_plane == 2) || (in_plane == 3));

            coord_data [0] [node] = double (x_plane);
            coord_data [1] [node] = double (y_plane);
            coord_data [2] [node] = double (z_plane);
        }

        Mesh_data::Coordinates<double>* coordinates = Mesh_data::Coordinates<double>::build_from (coord_data);

        // Side sets
        std::vector<Mesh_data::Side_set*> side_sets(0);

        // Node sets
        std::vector<Mesh_data::Node_set*> node_sets(0);

        // Data container
        data = Mesh_data::Data::build_from(parameters, coordinates, element_blocks, side_sets, node_sets);

    }

};

#endif /* _EXAMPLE_MESH_H_ */

