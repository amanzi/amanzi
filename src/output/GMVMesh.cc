/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "AmanziVector.hh"
#include "GMVMesh.hh"

namespace Amanzi {
namespace GMV {

static inline void
write_mesh_to_file_(const AmanziMesh::Mesh& mesh, std::string filename)
{
  int dim = mesh.getSpaceDimension();
  int dim_cell = mesh.getManifoldDimension();
  gmvwrite_openfile_ir_ascii((char*)filename.c_str(), 4, 8);

  // Write node info
  unsigned int num_nodes =
    mesh.getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
  double* x = new double[num_nodes];
  double* y = new double[num_nodes];
  double* z = new double[num_nodes];

  AmanziGeometry::Point xc(dim);
  for (int i = 0; i < num_nodes; i++) {
    xc = mesh.getNodeCoordinate(i);
    x[i] = xc[0];
    y[i] = xc[1];
    if (dim == 3)
      z[i] = xc[2];
    else
      z[i] = 0.0;
  }
  gmvwrite_node_data(&num_nodes, x, y, z);

  delete[] z;
  delete[] y;
  delete[] x;

  // Write cell info
  unsigned int num_cells =
    mesh.getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  gmvwrite_cell_header(&num_cells);

  int max_nodes = mesh.getCellMaxNodes();
  unsigned int* xh = new unsigned int[max_nodes];

  int max_faces = mesh.getCellMaxFaces();
  int* nverts = new int[max_faces];
  unsigned int* yh = new unsigned int[72];

  for (int i = 0; i < num_cells; i++) {
    auto cnodes = mesh.getCellNodes(i);

    int nnodes = cnodes.size();
    for (int j = 0; j < nnodes; j++) xh[j] = cnodes[j] + 1;

    if (dim == 3 && dim_cell == 3) {
      if (nnodes == 8) {
        gmvwrite_cell_type((char*)"phex8", 8, xh);
      } else {
        auto [cfaces, fdirs] = mesh.getCellFacesAndDirections(i);
        int nfaces = cfaces.size();

        for (int j = 0, n = 0; j < nfaces; j++) {
          auto fnodes = mesh.getFaceNodes(cfaces[j]);
          nverts[j] = fnodes.size();
          for (int k = 0; k < nverts[j]; k++) yh[n++] = fnodes[k] + 1;
        }
        gmvwrite_general_cell_type((char*)"general", nverts, nfaces, yh);
      }
    } else if (dim == 3 && dim_cell == 2) {
      gmvwrite_cell_type((char*)"general 1", nnodes, xh);
    } else if (dim == 2) {
      if (nnodes == 4)
        gmvwrite_cell_type((char*)"quad", 4, xh);
      else
        gmvwrite_cell_type((char*)"general 1", nnodes, xh);
    }
  }

  delete[] xh;
  delete[] nverts;
  delete[] yh;
}


void
create_mesh_file(const AmanziMesh::Mesh& mesh, std::string filename)
{
  write_mesh_to_file_(mesh, filename);
  gmvwrite_closefile();
}


void
open_data_file(std::string meshfile,
               std::string filename,
               unsigned int num_nodes,
               unsigned int num_cells)
{
  gmvwrite_openfile_ir_ascii((char*)filename.c_str(), 4, 8);
  gmvwrite_nodes_fromfile((char*)meshfile.c_str(), num_nodes);
  gmvwrite_cells_fromfile((char*)meshfile.c_str(), num_cells);
}


void
suffix_no(std::string& suffix, unsigned int cycleno)
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
    suffix[digits - i] = '0' + (cycleno / div) % 10;
    div *= 10;
  }
}


void
open_data_file(std::string meshfile,
               std::string filename,
               unsigned int num_nodes,
               unsigned int num_cells,
               unsigned int cycleno,
               unsigned int digits)
{
  std::string suffixstr(digits + 1, '.');
  suffix_no(suffixstr, cycleno);
  filename.append(suffixstr);

  gmvwrite_openfile_ir_ascii((char*)filename.c_str(), 4, 8);
  gmvwrite_nodes_fromfile((char*)meshfile.c_str(), num_nodes);
  gmvwrite_cells_fromfile((char*)meshfile.c_str(), num_cells);
}


void
open_data_file(const AmanziMesh::Mesh& mesh, const std::string& filename)
{
  write_mesh_to_file_(mesh, filename);
}


void
open_data_file(const AmanziMesh::Mesh& mesh,
               std::string filename,
               unsigned int cycleno,
               unsigned int digits)
{
  std::string suffixstr(digits + 1, '.');
  suffix_no(suffixstr, cycleno);
  filename.append(suffixstr);
  write_mesh_to_file_(mesh, filename);
}


void
write_time(const double time)
{
  gmvwrite_probtime(time);
}


void
write_cycle(const int cycle)
{
  gmvwrite_cycleno(cycle);
}


void
start_data()
{
  gmvwrite_variable_header();
}


void
write_node_data(const Vector_type& x, std::string varname)
{
  auto node_data = x.getData();
  gmvwrite_variable_name_data(1, (char*)varname.c_str(), (void*)node_data.get());
}


void
write_node_data(const MultiVector_type& x, const unsigned int component, std::string varname)
{
  auto node_data = x.getData(component);
  gmvwrite_variable_name_data(1, (char*)varname.c_str(), (void*)node_data.get());
}


void
write_cell_data(const Vector_type& x, std::string varname)
{
  auto cell_data = x.getData();
  gmvwrite_variable_name_data(0, (char*)varname.c_str(), (void*)cell_data.get());
}


void
write_cell_data(const MultiVector_type& x, const unsigned int component, std::string varname)
{
  auto cell_data = x.getData(component);
  gmvwrite_variable_name_data(0, (char*)varname.c_str(), (void*)cell_data.get());
}


void
write_face_data(const Vector_type& x, std::string varname)
{
  auto face_data = x.getData();
  gmvwrite_variable_name_data(2, (char*)varname.c_str(), (void*)face_data.get());
}


void
close_data_file()
{
  gmvwrite_variable_endvars();
  gmvwrite_closefile();
}

} // namespace GMV
} // namespace Amanzi
