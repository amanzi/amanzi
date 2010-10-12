#include "Side_set.hh"
#include "dbc.hh"

namespace Mesh_data
{

bool Side_set::valid () const
{

    bool result = (num_sides_ > 0);
    // result &= (num_nodes_ > 0);

    result &= (element_list_.size ()    == num_sides_);
    result &= (side_list_.size ()       == num_sides_);
    // result &= (node_list_.size ()       == num_nodes_);
    // result &= (node_count_list_.size () == num_sides_);
    // result &= (node_factors_.size ()    == num_nodes_ || node_factors_.size () == 0);

    return result;

}

void Side_set::take_data_from (std::vector<int>& element_list,
                               std::vector<int>& side_list,
                               std::string& name
                               // std::vector<int>& node_list,
                               // std::vector<int>& node_count_list,
                               // std::vector<double>& node_factors
                               )
{
    num_sides_ = element_list.size();
    // num_nodes_ = node_list.size();

    element_list_.swap (element_list);
    side_list_.swap (side_list);
    name_.swap (name);
    // node_list_.swap (node_list);
    // node_count_list_.swap (node_count_list);
    // node_factors_.swap (node_factors);

    ASSERT (valid ());

}

Side_set* Side_set::build_from (int id,
                                std::vector<int>& element_list,
                                std::vector<int>& side_list,
                                std::string& name
                                // std::vector<int>& node_list,
                                // std::vector<int>& node_count_list,
                                // std::vector<double>& node_factors
                                )
{

    Side_set* set = new Side_set (id);

    set->take_data_from (element_list, side_list, name // , node_list, node_count_list, node_factors
                         );

    return set;

}


void Side_set::to_stream (std::ostream& stream, bool verbose) const
{

    stream << "Side set: " << set_id_ << "\n";
    stream << "  Name: "   << name_   << "\n";
    stream << "  Number of sides: " << num_sides_ << "\n";
    // stream << "  Number of nodes: " << num_nodes_ << "\n";
    if (verbose)
    {
        stream << "  Elements and Side of Elements:\n";
        for (int element = 0; element < num_sides_ ; ++element)
            stream << "    " << element
                   << ".  "  << element_list_ [element]
                   << "  "   << side_list_ [element] << "\n";

    }


}


}
