#include "Element_block.hh"

#include <string.h>
#include "dbc.hh"

namespace Mesh_data
{

bool Element_block::valid () const
{
    bool valid = ( num_elements_ > 0);
    valid &= (num_nodes_per_element_ > 0);
    valid &= (num_attributes_ >= 0);
    valid &= (connectivity_map_.size () == num_elements_ * num_nodes_per_element_);
    valid &= (attribute_map_.size ()    == num_elements_ * num_attributes_);

    return valid;
}


std::vector<int> Element_block::connectivity (int element) const
{
    std::vector<int> nodes (num_nodes_per_element_);
    connectivity (element, nodes.begin ());
    return nodes;
}


void Element_block::take_data_from (int num_elements, 
                                    std::string name,
                                    ELEMENT_TYPE element_type,
                                    std::vector<int>& connectivity_map,
                                    std::vector<double>& attribute_map)
{
    element_type_ = element_type;
    name_ = name;

    num_elements_          = num_elements;
    num_nodes_per_element_ = connectivity_map.size () / num_elements;
    num_attributes_        = attribute_map.size ()    / num_elements;

    std::swap (connectivity_map_, connectivity_map);
    std::swap (attribute_map_,    attribute_map);

    ASSERT (valid ());
}


Element_block* Element_block::build_from (int id,
                                          std::string name,
                                          int num_elements,
                                          ELEMENT_TYPE element_type,
                                          std::vector<int>& connectivity_map,
                                          std::vector<double>& attribute_map)
{
    Element_block* element_block = new Element_block (id);
    element_block->take_data_from (num_elements, name, element_type, connectivity_map, attribute_map);
    return element_block;
}


void Element_block::to_stream (std::ostream& stream, bool verbose) const
{
    stream << "Element block " << block_id_ << ":\n";
    stream << "  Element type: " << type_to_name (element_type_) << "\n";
    stream << "  Number of elements: " <<  num_elements_ << "\n";
}




}
