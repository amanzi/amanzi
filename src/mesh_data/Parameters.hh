#ifndef _EXODUS_PARAMS_
#define _EXODUS_PARAMS_

#include <vector>
#include <string>

namespace Mesh_data
{

struct Parameters
{

    Parameters (std::string, int dimensions, int num_nodes, int num_elements,
                int num_element_blocks, int num_node_sets, int num_side_sets,
                const std::vector<int>& element_block_ids,
                const std::vector<int>& node_set_ids,
                const std::vector<int>& side_set_ids);

    std::string title_;
    unsigned int dimensions_;
    unsigned int num_nodes_;
    unsigned int num_elements_;

    unsigned int num_element_blocks_;
    unsigned int num_node_sets_;
    unsigned int num_side_sets_;

    std::vector<int> element_block_ids_;
    std::vector<int> node_set_ids_;
    std::vector<int> side_set_ids_;

    bool valid () const;
    
    int dimensions () const { return dimensions_; }

    void to_stream (std::ostream& stream) const;

    bool ok_node_id (int node)        const { return (node    >= 0) && (node    < num_nodes_);          }
    bool ok_element_id (int element)  const { return (element >= 0) && (element < num_elements_);       }
    bool ok_node_set (int set)        const { return (set     >= 0) && (set     < num_node_sets_);      }
    bool ok_side_set (int set)        const { return (set     >= 0) && (set     < num_side_sets_);      }
    bool ok_element_block (int block) const { return (block   >= 0) && (block   < num_element_blocks_); }


    virtual ~Parameters () { }


};

}

#endif
