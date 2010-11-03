#ifndef _GMV_MESH_
#define _GMV_MESH_

#include <string>
//#include "Mesh.hh"
#include "Mesh_maps_base.hh"
#include "Epetra_Vector.h"

extern "C" {
#include "gmvwrite.h"
}

namespace GMV {
  enum Data_type {
    CELL,
    NODE,
    FACE
  };

  // Write a GMV file containing only mesh data to be used as a "fromfile".
  void create_mesh_file(Mesh_maps_base &mesh_maps, std::string filename);

  // Opens and initializes a GMV file for writing which references a "fromfile" for mesh definition.
  void open_data_file(std::string mesh_fromfile, std::string filename, unsigned int num_nodes, unsigned int num_cells);
  // Opens and initializes a GMV file which contains mesh data, i.e. doesn't use a "fromfile".
  void open_data_file(Mesh_maps_base &mesh_maps, std::string filename);

  void start_variables();

  // Writes node data to files which has previously been opened with open_data_file.
  void write_node_data(const Epetra_Vector &x, std::string varname);

  // Writes cell data to files which has previously been opened with open_data_file.
  void write_cell_data(const Epetra_Vector &x, std::string varname);

  // Writes the cycle number and time
  void write_cycle (const int cycle);
  void write_time (const double time);

  // Finalizes a GMV file which has previously been opened with open_data_file.
  void close_data_file();
}

#endif  /* _GMV_MESH_ */
