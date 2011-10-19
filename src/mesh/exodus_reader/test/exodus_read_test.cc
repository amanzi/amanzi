#include <exodusII.h>
#include "dbc.hh"

#include <iostream>
#include <vector>
#include <iterator>

int main(int argc, char *argv[])
{
    
    int ret_val;
    int exodus_word_size (0);
    int system_word_size (sizeof (double));
    float version (0.0);

    ASSERT (argc > 1);

    int id = ex_open (argv [1], EX_READ, &system_word_size, &exodus_word_size, &version);

    std::cout << "Exodus id: " << id << std::endl;
    std::cout << "System word size: " << system_word_size << std::endl;
    std::cout << "Exodus file word size: " << exodus_word_size << std::endl;
    std::cout << "Version number: " << version << std::endl;

    ASSERT (id > 0);

    int dimension (0);
    int num_nodes (0);
    int num_elements (0);
    int num_element_blocks (0);
    int num_node_sets (0);
    int num_side_sets (0);

    char title_data [MAX_LINE_LENGTH];

    ret_val = ex_get_init (id, title_data, 
                           &dimension, &num_nodes, &num_elements, 
                           &num_element_blocks, &num_node_sets, &num_side_sets);
    ASSERT (ret_val >= 0);

    std::cout << "Dimension: "                << dimension          << std::endl;
    std::cout << "Title: "                    << title_data         << std::endl;
    std::cout << "Number of nodes: "          << num_nodes          << std::endl;
    std::cout << "Number of elements: "       << num_elements       << std::endl;
    std::cout << "Number of element blocks: " << num_element_blocks << std::endl;
    std::cout << "Number of node sets: "      << num_node_sets      << std::endl;
    std::cout << "Number of side sets: "      << num_side_sets      << std::endl;

    std::vector<int> side_set_ids (num_side_sets);
    ret_val = ex_get_side_set_ids (id, &side_set_ids [0]);
    ASSERT (ret_val >= 0);

    std::cout << "Side set ids: ";
    std::copy (side_set_ids.begin (), side_set_ids.end (), std::ostream_iterator<int>(std::cout, ", "));
    std::cout << std::endl;

    for (std::vector<int>::const_iterator side_set = side_set_ids.begin ();
         side_set != side_set_ids.end ();
         ++side_set)
    {
        std::cout << "Loading side set: " << *side_set << std::endl;

        int num_sides (0);
        int num_dist_factors (0);
        
        ret_val = ex_get_side_set_param (id, *side_set, &num_sides, &num_dist_factors);
        ASSERT (ret_val >= 0);

        std::cout << "  Number of sides: "                << num_sides        << std::endl;
        std::cout << "  Number of distribution factors: " << num_dist_factors << std::endl;

    }

    return 0;
    

}
