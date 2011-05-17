#ifndef _EXODUS_READERS_HH_
#define _EXODUS_READERS_HH_

#include "Data.hh"
#include "Exodus_file.hh"
#include "Exodus_error.hh"
#include "Coordinates.hh"
#include "Element_block.hh"
#include "Side_set.hh"
#include "Node_set.hh"

namespace Amanzi {
namespace Exodus {

AmanziMesh::Data::Parameters* read_parameters (Exodus_file file);

AmanziMesh::Data::Coordinates<double>* read_coordinates (Exodus_file file, int num_nodes, int dimensions);

AmanziMesh::Data::Element_block* read_element_block (Exodus_file file, int block_id);

AmanziMesh::Data::Side_set* read_side_set (Exodus_file file, int set_id);

AmanziMesh::Data::Node_set* read_node_set (Exodus_file file, int set_id);

AmanziMesh::Data::Data* read_exodus_file (const char *);

AmanziMesh::Data::Data* read_exodus_file (const Exodus_file& file);

AmanziMesh::Cell_type read_element_type (const char *);

} // namespace Exodus
} // namespace Amanzi


#endif
