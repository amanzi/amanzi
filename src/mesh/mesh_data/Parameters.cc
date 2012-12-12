#include "Parameters.hh"

#include <iterator>

#include "Utils.hh"
#include "dbc.hh"

namespace Amanzi {
namespace AmanziMesh {
namespace Data {

bool Parameters::valid () const
{

    bool result = (dimensions_ > 0);
    result &= (dimensions_ <= 3);
    result &= (num_nodes_ > 0);
    result &= (num_elements_ > 0);
    result &= (num_element_blocks_ > 0);

    result &= (element_block_ids_.size () == num_element_blocks_);
    result &= (node_set_ids_.size () == num_node_sets_);
    result &= (side_set_ids_.size () == num_side_sets_);

    return result;

}


Parameters::Parameters (std::string title,
                        int dimensions, int num_nodes, int num_elements, 
                        int num_element_blocks, int num_node_sets, int num_side_sets,
                        const std::vector<int>& element_block_ids,
                        const std::vector<int>& node_set_ids,
                        const std::vector<int>& side_set_ids) : 
    title_ (title),
    dimensions_ (dimensions),
    num_nodes_ (num_nodes),
    num_elements_ (num_elements),
    num_element_blocks_ (num_element_blocks),
    num_node_sets_ (num_node_sets),
    num_side_sets_ (num_side_sets),
    element_block_ids_ (element_block_ids),
    node_set_ids_ (node_set_ids),
    side_set_ids_ (side_set_ids)
    
{
    ASSERT (valid ());
}


void Parameters::to_stream (std::ostream& stream, const bool& verbose) const
{
    stream << "Mesh Parameters:\n";
    stream << "  Title: " << title_ << "\n";
    stream << "  Dimension: " << dimensions_ << "\n";
    stream << "  Number of Nodes: " << num_nodes_ << "\n";
    stream << "  Number of Elements: " << num_elements_ << "\n";
    stream << "  Number of Element Blocks: " << num_element_blocks_ << "\n";
    stream << "  Number of Node Sets: " << num_node_sets_ << "\n";
    stream << "  Number of Side Sets: " << num_side_sets_ << "\n";

    if (num_element_blocks_ > 0)
    {
        stream << "  Element Block Ids: ";
        std::copy (element_block_ids_.begin (), element_block_ids_.end (), std::ostream_iterator<int>(stream, ", "));
        stream << "\n";
    }

    if (num_node_sets_ > 0)
    {
        stream << "  Node Set Ids: ";
        std::copy (node_set_ids_.begin (), node_set_ids_.end (), std::ostream_iterator<int>(stream, ", "));
        stream << "\n";
    }

    if (num_side_sets_ > 0)
    {
        stream << "  Side Set Ids: ";
        std::copy (side_set_ids_.begin (), side_set_ids_.end (), std::ostream_iterator<int>(stream, ", "));
        stream << "\n";
    }

}


} // namespace Data
} // namespace Mesh
} // namespace Amanzi
