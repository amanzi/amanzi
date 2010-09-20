#ifndef _ELEMENT_BLOCK_HH_
#define _ELEMENT_BLOCK_HH_

#include <iostream>
#include <vector>
#include <string>

#include "Element_types.hh"
#include "dbc.hh"

namespace Mesh_data
{


class Element_block
{

private:

    int block_id_;
    std::string name_;

    int num_elements_;
    int num_nodes_per_element_;
    int num_attributes_;

    ELEMENT_TYPE element_type_;
    std::vector<int> connectivity_map_;
    std::vector<double> attribute_map_;

    Element_block (int block_id, std::string name = "") : block_id_ (block_id), name_ (name) { }

    bool valid () const;

public:

    // Block information
    // -----------------

    ELEMENT_TYPE element_type () const  { return element_type_; }
    int nodes_per_element () const      { return num_nodes_per_element_; }
    int num_elements () const           { return num_elements_; }
    bool has_attributes () const        { return num_attributes_ > 0; }
    bool ok_element (int element) const { return (element >= 0) && (element < num_elements_); }
    const std::string& name () const    { return name_; }

    // Connectivity information
    // ------------------------

    const std::vector<int>& connectivity () const { return connectivity_map_; }
    template <typename IT> void connectivity (int element, IT storage) const;
    std::vector<int> connectivity (int element) const;

    void to_stream (std::ostream& stream, bool verbose = false) const;


    // Construction and Data assignment.
    // --------------------------------

    void take_data_from (int num_elements, 
                         std::string name,
                         ELEMENT_TYPE type, 
                         std::vector<int>& connectivity_map, 
                         std::vector<double>& attribute_map);
    
    static Element_block* build_from (int id, 
                                      std::string name,
                                      int num_elements, 
                                      ELEMENT_TYPE type, 
                                      std::vector<int>& connectivity_map, 
                                      std::vector<double>& attribute_map);
    
    virtual ~Element_block () { };

};

template <typename IT>
void Element_block::connectivity (int element, IT storage) const
{
    ASSERT (ok_element (element));

    const std::vector<int>::const_iterator start = 
        connectivity_map_.begin () + element * num_nodes_per_element_;

    std::copy (start, start + num_nodes_per_element_, storage);

}

}
#endif
