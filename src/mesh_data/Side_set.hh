#ifndef _SIDE_SET_HH_
#define _SIDE_SET_HH_

#include <vector>
#include <iostream>

namespace Mesh_data
{

class Side_set
{

    int num_sides_;
    int num_nodes_;

    int set_id_;

    std::vector<int> element_list_;
    std::vector<int> side_list_;

    std::vector<int> node_list_;
    std::vector<int> node_count_list_;
    std::vector<double> node_factors_;

    bool valid () const;

public:

    Side_set (int set_id) : set_id_ (set_id) { }

    bool has_node_factors () const { return node_factors_.size () > 0; }

    int num_sides () const { return num_sides_; }
    int num_nodes () const { return num_nodes_; }

    const std::vector<int> element_list ()    const { return element_list_;    }
    const std::vector<int> side_list ()       const { return side_list_;       }
    const std::vector<int> node_list ()       const { return node_list_;       }
    const std::vector<int> node_count_list () const { return node_count_list_; }
    const std::vector<double> node_factors () const { return node_factors_;    }

    int id () const { return set_id_; }

    void to_stream (std::ostream& stream, bool verbose = false) const;

    void take_data_from (std::vector<int>& element_list,
                         std::vector<int>& side_list,
                         std::vector<int>& node_list,
                         std::vector<int>& node_count_list,
                         std::vector<double>& node_factors);

    static Side_set* build_from (int id, 
                                 std::vector<int>& element_list,
                                 std::vector<int>& side_list,
                                 std::vector<int>& node_list,
                                 std::vector<int>& node_count_list,
                                 std::vector<double>& node_factors);



};

}

#endif
