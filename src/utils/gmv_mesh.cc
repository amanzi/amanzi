#include "gmv_mesh.hh"

namespace GMV {

  static inline void write_mesh_to_file_(Mesh_maps_base &mesh_map, std::string filename)
  {
    gmvwrite_openfile_ir_ascii((char*)filename.c_str(), 4, 8);

    // Write node info
    unsigned int num_nodes = mesh_map.count_entities(Mesh_data::NODE, OWNED);
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
  
    delete [] z;
    delete [] y;
    delete [] x;

    // Write cell info
    unsigned int num_cells = mesh_map.count_entities(Mesh_data::CELL, OWNED);

    gmvwrite_cell_header(&num_cells);

    unsigned int *xh = new unsigned int[8];
    for (int i=0; i<num_cells; i++) {
      mesh_map.cell_to_nodes(i,xh,xh+8);
      for (int j=0; j<8; j++) xh[j]++;
      gmvwrite_cell_type((char*) "phex8",8,xh);
    }

    delete [] xh;
    
  }

  void create_mesh_file(Mesh_maps_base &mesh_map, std::string filename)
  {
    write_mesh_to_file_(mesh_map, filename);
    gmvwrite_closefile();
  }

  void open_data_file(std::string meshfile, std::string filename, unsigned int num_nodes, unsigned int num_cells) {

    gmvwrite_openfile_ir_ascii((char*)filename.c_str(), 4, 8);
    gmvwrite_nodes_fromfile((char*) meshfile.c_str(), num_nodes);
    gmvwrite_cells_fromfile((char*) meshfile.c_str(), num_cells);
  }


  void suffix_no(std::string &suffix, unsigned int cycleno) {
    
    unsigned int digits = suffix.length() - 1;

    suffix[0] = '.';

    if (cycleno >= pow(10.0, static_cast<double>(digits)))  throw std::exception();

    // suffix[digits] = '0' + cycleno%10;
    
    // int div = 10;
    // for (int i=1; i<digits; i++) {
    //   suffix[digits-i] = '0' + (cycleno/div)%div;
    //   div *=10;
    // }  

    int div = 1;
    for(int i = 0; i < digits; ++i) {
      suffix[digits-i] = '0' + (cycleno/div)%10; 
      div *= 10;
    }

  }
  


  void open_data_file(std::string meshfile, std::string filename, unsigned int num_nodes, unsigned int num_cells, unsigned int cycleno, unsigned int digits) {
    
    string suffixstr(digits+1,'.');
    
    suffix_no(suffixstr, cycleno);
    filename.append(suffixstr);
    
    gmvwrite_openfile_ir_ascii((char*)filename.c_str(), 4, 8);
    gmvwrite_nodes_fromfile((char*) meshfile.c_str(), num_nodes);
    gmvwrite_cells_fromfile((char*) meshfile.c_str(), num_cells);
  }

  void open_data_file(Mesh_maps_base &mesh_map, std::string filename) {

    unsigned int num_nodes = mesh_map.count_entities(Mesh_data::NODE, OWNED);
    unsigned int num_cells = mesh_map.count_entities(Mesh_data::CELL, OWNED);

    write_mesh_to_file_(mesh_map, filename);
  }

  void open_data_file(Mesh_maps_base &mesh_map, std::string filename, unsigned int cycleno, unsigned int digits) {

    unsigned int num_nodes = mesh_map.count_entities(Mesh_data::NODE, OWNED);
    unsigned int num_cells = mesh_map.count_entities(Mesh_data::CELL, OWNED);

    string suffixstr(digits+1,'.');
    
    suffix_no(suffixstr, cycleno);
    filename.append(suffixstr); 
    
    write_mesh_to_file_(mesh_map, filename);
  }
  

  void write_time(const double time) {
    gmvwrite_probtime(time);
  }

  void write_cycle(const int cycle) {
    gmvwrite_cycleno(cycle);
  }

  void start_data() {
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
