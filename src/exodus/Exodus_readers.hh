#ifndef _EXODUS_READERS_HH_
#define _EXODUS_READERS_HH_

#include "Data.hh"
#include "Exodus_file.hh"
#include "Exodus_error.hh"
#include "Coordinates.hh"
#include "Element_block.hh"
#include "Side_set.hh"
#include "Node_set.hh"

namespace ExodusII
{

Mesh_data::Parameters* read_parameters (Exodus_file file);

Mesh_data::Coordinates<double>* read_coordinates (Exodus_file file, int num_nodes, int dimensions);

Mesh_data::Element_block* read_element_block (Exodus_file file, int block_id);

Mesh_data::Side_set* read_side_set (Exodus_file file, int set_id);

Mesh_data::Node_set* read_node_set (Exodus_file file, int set_id);

Mesh_data::Data* read_exodus_file (const char *);

Mesh_data::Data* read_exodus_file (const Exodus_file& file);

Mesh_data::ELEMENT_TYPE read_element_type (const char *);

}



#endif
