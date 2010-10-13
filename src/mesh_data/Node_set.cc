#include "Node_set.hh"

#include "dbc.hh"

namespace Mesh_data
{

void Node_set::take_data_from (std::vector<int>& node_list, std::vector<double>& dist_factors, std::string& name)
{
    num_nodes_ = node_list.size ();
    num_dist_factors_ = dist_factors.size ();

    node_list_.swap (node_list);
    node_dist_factors_.swap (dist_factors);

    name_.swap (name);

    ASSERT (valid ());
}


Node_set* Node_set::build_from (int set_id, std::vector<int>& node_list, std::vector<double>& dist_factors,
                                std::string& name)
{
    Node_set *set = new Node_set (set_id);
    set->take_data_from (node_list, dist_factors, name);
    return set;
}

bool Node_set::valid () const
{
    return (node_list_.size () == num_nodes_) && (node_dist_factors_.size () == num_dist_factors_);
}

void Node_set::to_stream (std::ostream& stream, bool verbose) const
{
    stream << "Node set: " << set_id_ << std::endl;
    stream << "  Number of nodes: " << num_nodes_ << std::endl;
    stream << "  Distribution factors: " << num_dist_factors_ << std::endl;
}

}
