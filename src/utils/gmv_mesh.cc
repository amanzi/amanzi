#include "gmv_mesh.hh"
#include "Entity_kind.hh"
#include "Element_category.hh"

namespace GMV {

  static inline void write_mesh_to_file_(STK_mesh::Mesh_maps_stk &mesh_map, std::string filename)
  {
    gmvwrite_openfile_ir_ascii((char*)filename.c_str(), 4, 8);

    // Write node info
    unsigned int num_nodes = mesh_map.count_entities(Mesh_data::NODE, STK_mesh::OWNED);
    double *x = new double [num_nodes];
    double *y = new double [num_nodes];
    double *z = new double [num_nodes];
    double xc[3];

    for (int i=0; i<num_nodes; i++) {
      mesh_map.node_to_coordinates(i,xc,xc+3);
      x[i] = xc[0];
      y[i] = xc[1];
      z[i] = xc[2];
    }
    gmvwrite_node_data(&num_nodes, x, y, z);
  
    delete x;
    delete y;
    delete z;

    // Write cell info
    unsigned int num_cells = mesh_map.count_entities(Mesh_data::CELL, STK_mesh::OWNED);

    gmvwrite_cell_header(&num_cells);

    int *xh = new int[8];
    for (int i=0; i<num_cells; i++) {
      mesh_map.cell_to_nodes(i,xh,xh+8);
      for (int j=0; j<8; j++) xh[j]++;
      gmvwrite_cell_type((char*) "phex8",8,xh);
    }
    
  }

  void create_mesh_file(STK_mesh::Mesh_maps_stk &mesh_map, std::string filename)
  {
    write_mesh_to_file_(mesh_map, filename);
    gmvwrite_closefile();
  }

  void open_data_file(std::string meshfile, std::string filename, unsigned int num_nodes, unsigned int num_cells) {

    gmvwrite_openfile_ir_ascii((char*)filename.c_str(), 4, 8);
    gmvwrite_nodes_fromfile((char*) meshfile.c_str(), num_nodes);
    gmvwrite_cells_fromfile((char*) meshfile.c_str(), num_cells);
    gmvwrite_variable_header();
  }

  void open_data_file(STK_mesh::Mesh_maps_stk &mesh_map, std::string filename) {

    unsigned int num_nodes = mesh_map.count_entities(Mesh_data::NODE, STK_mesh::OWNED);
    unsigned int num_cells = mesh_map.count_entities(Mesh_data::CELL, STK_mesh::OWNED);

    write_mesh_to_file_(mesh_map, filename);
    gmvwrite_variable_header();
  }

  void write_node_data(const Epetra_Vector &x, std::string varname) {
    double *node_data;
    int err = x.ExtractView(&node_data);
    gmvwrite_variable_name_data(NODE, (char *) varname.c_str(), node_data);
  }

  void write_cell_data(const Epetra_Vector &x, std::string varname) {
    double *cell_data;
    int err = x.ExtractView(&cell_data);
    gmvwrite_variable_name_data(CELL, (char *) varname.c_str(), cell_data);
  }

  void close_data_file() {
    gmvwrite_variable_endvars();
    gmvwrite_closefile();
  }
}
