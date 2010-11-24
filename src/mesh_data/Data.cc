#include "Data.hh"

#include "Utils.hh"

namespace Mesh_data
{

void Data::take_data_from (Parameters* params,
                           Coordinates<double>* coords,
                           std::vector<Element_block*> blocks,
                           std::vector<Side_set*> side_sets,
                           std::vector<Node_set*> node_sets)
{

    coords_ = std::auto_ptr<Coordinates<double> >(coords);
    params_ = std::auto_ptr<Parameters>(params);
    element_blocks_.swap (blocks);
    side_sets_.swap (side_sets);
    node_sets_.swap (node_sets);

}

Data* Data::build_from (Parameters* params,
                        Coordinates<double>* coords,
                        std::vector<Element_block*> blocks,
                        std::vector<Side_set*> side_sets,
                        std::vector<Node_set*> node_sets)
{

    Data* mesh = new Data ();

    mesh->take_data_from (params, coords, blocks, side_sets, node_sets);

    return mesh;

}

Data::~Data ()
{

    Utils::reclaim (element_blocks_);
    Utils::reclaim (side_sets_);
    Utils::reclaim (node_sets_);

}

void Data::to_stream (std::ostream& stream, const bool& verbose) const
{
    params_->to_stream (stream, verbose);

    for (int element = 0; element < element_blocks (); ++element)
        element_block (element).to_stream (stream, verbose);

    for (int set = 0; set < side_sets (); ++set)
        side_set (set).to_stream (stream, verbose);

    for (int set = 0; set < node_sets (); ++set)
        node_set (set).to_stream (stream, verbose);

}

const Element_block& Data::element_block (int id) const
{
    ASSERT (params_->ok_element_block (id));
    return *element_blocks_ [id];
}

const Side_set& Data::side_set (int id) const
{
    ASSERT (params_->ok_side_set (id));
    return *side_sets_ [id];
}

const Node_set& Data::node_set (int id) const
{
    ASSERT (params_->ok_node_set (id));
    return *node_sets_ [id];
}

}
