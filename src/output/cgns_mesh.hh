#ifndef _CGNS_MESH_
#define _CGNS_MESH_

#include <string>
#include <vector>
#include "Mesh_maps_base.hh"
#include "Epetra_Vector.h"

extern "C" {
#include "cgnslib.h"
}


namespace CGNS {
    
  // Write a CGNS file containing only mesh data to be used as a "fromfile".
  void create_mesh_file(Mesh_maps_base &mesh_maps, std::string filename);

  // Opens and initializes a CGNS file for writing which references a "fromfile" for mesh definition.
  void open_data_file(std::string mesh_fromfile, std::string filename,  
		      unsigned int num_nodes, unsigned int num_cells);
  // Opens an CGNS file which already contains mesh data or link to mesh "fromfile".
  void open_data_file(std::string filename);

  // Create timestep node in CGNS structure of which has previously been opened with open_data_file.
  void create_timestep(const double time, const int iter, Mesh_data::Entity_kind kind);
    
  // Writes data to files which has previously been opened with open_data_file.
  void write_field_data(const Epetra_Vector &x, std::string varname);
    
  // Writes data to files which has previously been opened with open_data_file.
    //void write_probe_data(const Epetra_Vector &x, std::string varname);
        
  // Finalizes a CGNS file which has previously been opened with open_data_file.
  void close_data_file();
}

#endif  /* _CGNS_MESH_ */
