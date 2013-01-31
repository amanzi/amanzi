#include "GMVMesh.hh"

namespace Amanzi {
namespace GMV {

static inline void write_mesh_to_file_(AmanziMesh::Mesh &mesh_map, std::string filename)
{
  int dim = mesh_map.space_dimension();
  gmvwrite_openfile_ir_ascii((char*)filename.c_str(), 4, 8);

  // Write node info
  unsigned int num_nodes = mesh_map.num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);
  double *x = new double[num_nodes];
  double *y = new double[num_nodes];
  double *z = new double[num_nodes];

  AmanziGeometry::Point xc(dim);
  for (int i=0; i<num_nodes; i++) {
    mesh_map.node_get_coordinates(i, &xc);
    x[i] = xc[0];
    y[i] = xc[1];
    if (dim == 3)  
      z[i] = xc[2];
    else 
      z[i] = 0.0;
  }
  gmvwrite_node_data(&num_nodes, x, y, z);
  
  delete [] z;
  delete [] y;
  delete [] x;

  // Write cell info
  unsigned int num_cells = mesh_map.num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  gmvwrite_cell_header(&num_cells);

  unsigned int *xh = new unsigned int[8];
  for (int i=0; i<num_cells; i++) {
    AmanziMesh::Entity_ID_List cnodes;
    mesh_map.cell_get_nodes(i, &cnodes);

    int nnodes = cnodes.size();
    for (int j=0; j<nnodes; j++) xh[j] = cnodes[j] + 1;

    if (dim == 3) {
      gmvwrite_cell_type((char*) "phex8", 8, xh);
    } else if (dim == 2) {
      if (nnodes == 4) 
        gmvwrite_cell_type((char*) "quad", 4, xh);
      else
        gmvwrite_cell_type((char*) "general 1", nnodes, xh);
    }
  }

  delete [] xh;
}


void create_mesh_file(AmanziMesh::Mesh &mesh_map, std::string filename)
{
  write_mesh_to_file_(mesh_map, filename);
  gmvwrite_closefile();
}


void open_data_file(std::string meshfile, std::string filename, unsigned int num_nodes, unsigned int num_cells) 
{
  gmvwrite_openfile_ir_ascii((char*) filename.c_str(), 4, 8);
  gmvwrite_nodes_fromfile((char*) meshfile.c_str(), num_nodes);
  gmvwrite_cells_fromfile((char*) meshfile.c_str(), num_cells);
}


void suffix_no(std::string &suffix, unsigned int cycleno) 
{
  unsigned int digits = suffix.length() - 1;

  suffix[0] = '.';
  if (cycleno >= pow(10.0, static_cast<double>(digits))) throw std::exception();

  // suffix[digits] = '0' + cycleno%10;
  // int div = 10;
  // for (int i=1; i<digits; i++) {
  //   suffix[digits-i] = '0' + (cycleno/div)%div;
  //   div *=10;
  // }  

  int div = 1;
  for (int i = 0; i < digits; ++i) {
    suffix[digits-i] = '0' + (cycleno/div) % 10; 
    div *= 10;
  }
}
  

void open_data_file(std::string meshfile, 
                    std::string filename, 
                    unsigned int num_nodes, 
                    unsigned int num_cells, 
                    unsigned int cycleno, 
                    unsigned int digits) 
{
  string suffixstr(digits+1,'.');
  suffix_no(suffixstr, cycleno);
  filename.append(suffixstr);
    
  gmvwrite_openfile_ir_ascii((char*) filename.c_str(), 4, 8);
  gmvwrite_nodes_fromfile((char*) meshfile.c_str(), num_nodes);
  gmvwrite_cells_fromfile((char*) meshfile.c_str(), num_cells);
}


void open_data_file(AmanziMesh::Mesh &mesh_map, std::string filename) 
{
  unsigned int num_nodes = mesh_map.num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);
  unsigned int num_cells = mesh_map.num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  write_mesh_to_file_(mesh_map, filename);
}


void open_data_file(AmanziMesh::Mesh &mesh_map, std::string filename, unsigned int cycleno, unsigned int digits) 
{
  unsigned int num_nodes = mesh_map.num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);
  unsigned int num_cells = mesh_map.num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  string suffixstr(digits+1,'.');
  suffix_no(suffixstr, cycleno);
  filename.append(suffixstr); 
  write_mesh_to_file_(mesh_map, filename);
}
  

void write_time(const double time) 
{
  gmvwrite_probtime(time);
}


void write_cycle(const int cycle) 
{
  gmvwrite_cycleno(cycle);
}


void start_data() 
{
  gmvwrite_variable_header();
}
    

void write_node_data(const Epetra_Vector &x, std::string varname) 
{
  double *node_data;
  int err = x.ExtractView(&node_data);
  gmvwrite_variable_name_data(AmanziMesh::NODE, (char*) varname.c_str(), node_data);
}


void write_cell_data(const Epetra_Vector &x, std::string varname) 
{
  double *cell_data;
  int err = x.ExtractView(&cell_data);
  gmvwrite_variable_name_data(0, (char*) varname.c_str(), cell_data);
}


void write_cell_data(const Epetra_MultiVector &x, const unsigned int component, std::string varname) 
{
  double **cell_data;
  int err = x.ExtractView(&cell_data);

  double *component_data = cell_data[component];
  gmvwrite_variable_name_data(0, (char*) varname.c_str(), component_data);
}


void write_face_data(const Epetra_Vector &x, std::string varname) 
{
  double *face_data;
  int err = x.ExtractView(&face_data);
  gmvwrite_variable_name_data(AmanziMesh::FACE, (char*) varname.c_str(), face_data);
}

    
void close_data_file() 
{
  gmvwrite_variable_endvars();
  gmvwrite_closefile();
}

}  // namespace GMV
}  // namespace Amanzi
