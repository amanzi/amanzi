#ifndef _NODE_SET_HH_
#define _NODE_SET_HH_

#include <vector>
#include <iostream>

namespace Mesh_data
{

class Node_set
{

    int set_id_;

    int num_nodes_;
    int num_dist_factors_;

    std::vector<int> node_list_;
    std::vector<double> node_dist_factors_;

    bool valid () const;

public:

    Node_set (int set_id) : set_id_ (set_id) { }

    bool has_node_factors () const { return num_dist_factors_ > 0; }

    int num_nodes () const { return num_nodes_; }

    int id () const { return set_id_; }

    const std::vector<int> node_list () const { return node_list_; }

    void to_stream (std::ostream& stream, bool verbose = false) const;

    void take_data_from (std::vector<int>& node_list, std::vector<double>& dist_factors);
    static Node_set* build_from (int set_id, std::vector<int>& node_list, std::vector<double>& dist_factors);


};


}


#endif
