#include "Exodus_readers.hh"
#include "Exodus_error.hh"

#include "Element_types.hh"

#include <exodusII.h>

namespace ExodusII
{

Mesh_data::Parameters* read_parameters (Exodus_file file)
{

    char title_data [MAX_LINE_LENGTH];

    // Global mesh information
    // -----------------------
    int dimensionality;
    int num_nodes;
    int num_elements;
    int num_element_blocks;
    int num_node_sets;
    int num_side_sets;
    int ret_val = ex_get_init (file.id, title_data, 
                               &dimensionality, &num_nodes, &num_elements, 
                               &num_element_blocks, &num_node_sets, &num_side_sets);

    if (ret_val < 0) throw ExodusII::ExodusError (ret_val);
    std::string title(title_data);

    std::vector<int> element_block_ids;
    std::vector<int> node_set_ids;
    std::vector<int> side_set_ids;

    // Element blocks
    // --------------
    if (num_element_blocks > 0)
    {
        element_block_ids.resize (num_element_blocks);
        ret_val = ex_get_elem_blk_ids (file.id, &element_block_ids [0]);
        if (ret_val < 0) throw ExodusII::ExodusError (ret_val);
    }


    // Node sets
    // ---------
    if (num_node_sets > 0)
    {
        node_set_ids.resize (num_node_sets);
        ret_val = ex_get_node_set_ids (file.id, &node_set_ids [0]);
        if (ret_val < 0) throw ExodusII::ExodusError (ret_val);
    }


    // Side sets
    // ---------
    if (num_side_sets > 0)
    {
        side_set_ids.resize (num_side_sets);
        ret_val = ex_get_side_set_ids (file.id, &side_set_ids [0]);
        if (ret_val < 0) throw ExodusII::ExodusError (ret_val);
    }



    
    Mesh_data::Parameters *params = new Mesh_data::Parameters (title, dimensionality, num_nodes, num_elements,
                                                     num_element_blocks, num_node_sets, num_side_sets,
                                                     element_block_ids, node_set_ids, side_set_ids);

    return params;

}


Mesh_data::Coordinates<double>* read_coordinates (Exodus_file file, int num_nodes, int dimensions)
{

    if (file.exodus_word_size == sizeof (float))
        return read_coordinates_impl_<float> (file, num_nodes, dimensions);

    if (file.exodus_word_size == sizeof (double))
        return read_coordinates_impl_<double> (file, num_nodes, dimensions);

    return NULL;
}


Mesh_data::Node_set* read_node_set (Exodus_file file, int set_id)
{

    int num_nodes;
    int num_dist_factors;
    int ret_val = ex_get_node_set_param (file.id, set_id, &num_nodes, &num_dist_factors);
    if (ret_val < 0) throw ExodusII::ExodusError (ret_val);

    std::vector<int> node_list(num_nodes);
    std::vector<double> node_dist_factors(num_dist_factors);

    ret_val = ex_get_node_set (file.id, set_id, &node_list [0]);
    if (ret_val < 0) throw ExodusII::ExodusError (ret_val);

    if (num_dist_factors > 0)
    {
        ret_val = ex_get_node_set_dist_fact (file.id, set_id, &node_dist_factors [0]);
        if (ret_val < 0) throw ExodusII::ExodusError (ret_val);
    }

    return Mesh_data::Node_set::build_from (set_id, node_list, node_dist_factors);

}


Mesh_data::Side_set* read_side_set (Exodus_file file, int set_id)
{

    int num_sides;
    int num_nodes;

    int ret_val = ex_get_side_set_param (file.id, set_id, &num_sides, &num_nodes);
    if (ret_val < 0) throw ExodusII::ExodusError (ret_val);


    std::vector<int> element_list(num_sides);
    std::vector<int> side_list(num_sides);

    std::vector<int> node_list(num_nodes);
    std::vector<int> node_count_list(num_sides);
    std::vector<double> node_factors(num_nodes);
        
    ret_val = ex_get_side_set (file.id, set_id, &element_list [0], &side_list [0]);
    if (ret_val < 0) throw ExodusII::ExodusError (ret_val);

    ret_val = ex_get_side_set_node_list (file.id, set_id, &node_count_list [0], &node_list [0]);
    if (ret_val < 0) throw ExodusII::ExodusError (ret_val);

    ret_val = ex_get_side_set_dist_fact (file.id, set_id, &node_factors [0]);
    if (ret_val < 0) throw ExodusII::ExodusError (ret_val);
    if (ret_val > 0)
    {
        std::vector<double>().swap (node_factors);
    }    

    return Mesh_data::Side_set::build_from (set_id, 
                                            element_list, 
                                            side_list, 
                                            node_list, 
                                            node_count_list, 
                                            node_factors);
    
}


Mesh_data::Element_block* read_element_block(Exodus_file file, int block_id)
{

    char element_name_data [MAX_STR_LENGTH];

    // Block size and type data
    int num_elements;
    int num_nodes_per_element;
    int num_attributes;
    int ret_val = ex_get_elem_block (file.id, block_id, element_name_data,
                                     &num_elements, &num_nodes_per_element, &num_attributes);
    if (ret_val < 0) throw ExodusII::ExodusError (ret_val);


    // Connectivity Map
    std::vector<int> connectivity_map(num_elements * num_nodes_per_element);
    ret_val = ex_get_elem_conn (file.id, block_id, &connectivity_map [0]);
    if (ret_val < 0) throw ExodusII::ExodusError (ret_val);

    // Attribute Map
    std::vector<double> attribute_map;
    if (num_attributes > 0)
    {
        attribute_map.resize (num_elements * num_attributes);
        ret_val = ex_get_elem_attr (file.id, block_id, &attribute_map [0]);
        if (ret_val < 0) throw ExodusII::ExodusError (ret_val);
    }

    Mesh_data::ELEMENT_TYPE element_type = read_element_type(element_name_data);

    return  Mesh_data::Element_block::build_from (block_id,
                                                  std::string (element_name_data),
                                                  num_elements,
                                                  element_type, 
                                                  connectivity_map,
                                                  attribute_map);

}

Mesh_data::Data* read_exodus_file (const char * filename)
{
    
    Exodus_file file (filename);

    Mesh_data::Parameters* parameters = read_parameters (file);

    const int num_nodes = parameters->num_nodes_;
    const int dimensions = parameters->dimensions_;

    Mesh_data::Coordinates<double>* coords = read_coordinates (file, num_nodes, dimensions);

    // Element blocks:
    std::vector<Mesh_data::Element_block*> element_blocks (parameters->num_element_blocks_);
    for (int i = 0; i < parameters->num_element_blocks_; ++i)
        element_blocks [i] = read_element_block (file, parameters->element_block_ids_ [i]);

    // Side sets:
    std::vector<Mesh_data::Side_set*> side_sets (parameters->num_side_sets_);
    for (int i = 0; i < parameters->num_side_sets_; ++i)
        side_sets [i] = read_side_set (file, parameters->side_set_ids_ [i]);

    // Node sets:
    std::vector<Mesh_data::Node_set*> node_sets (parameters->num_node_sets_);
    for (int i = 0; i < parameters->num_node_sets_; ++i)
        node_sets [i] = read_node_set (file, parameters->node_set_ids_ [i]);

    
    Mesh_data::Data* mesh = Mesh_data::Data::build_from (parameters, coords, element_blocks, side_sets, node_sets);

    return mesh;
        

}

Mesh_data::ELEMENT_TYPE read_element_type (const char * name)
{

    std::string id (name, name+3);

    if (id == "CIR")
        return Mesh_data::CIRCLE;
    if (id == "SPH")
        return Mesh_data::SPHERE;
    if (id == "TRU")
        return Mesh_data::TRUSS;
    if (id == "BEA")
        return Mesh_data::BEAM;
    if (id == "TRI")
        return Mesh_data::TRIANGLE;
    if (id == "QUA")
        return Mesh_data::QUAD;
    if (id == "SHE")
        return Mesh_data::SHELL;
    if (id == "TET")
        return Mesh_data::TETRA;
    if (id == "WED")
        return Mesh_data::WEDGE;
    if (id == "HEX")
        return Mesh_data::HEX;

    return Mesh_data::UNKNOWN;

}


}
