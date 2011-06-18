#include "cgns_mesh.hh"
//#include <vector>



// Handles to carry around: - keep here or in cgnswrite??
//struct cgns_fh
int file_idx;    // file index 
int base_idx;    // root node
int zone_idx;    // mesh zone
int coord_idx;   // node coordinates
int soln_idx;    // solution node index
int field_idx;    // solution node index
int section_idx;  // cell connectivity
int num_soln = 0;   // number of solution nodes
std::vector<double> timeList;
std::vector<string> nameList;

namespace Amanzi {
namespace CGNS {
    
void 
create_mesh_file(AmanziMesh::Mesh &mesh_maps, std::string filename)
{
  int isize[3][1];
  int nelem_start, nelem_end, nbdyelem;
	
	
  // open file
  //cg_set_file_type(CG_FILE_ADF);
  cg_open ((const char*)filename.c_str(), CG_MODE_WRITE, &file_idx);
	
  // create base
  int icelldim = 3;  // cell dimension: 3-volume, 2-face, ...
  int iphysdim = 3;  // physical dimension: 1-1D, 2-2D, 3-3D
  cg_base_write(file_idx, "Base", icelldim, iphysdim, &base_idx);
			
  // get num_nodes, num_cells
  unsigned int num_nodes = mesh_maps.num_entities(AmanziMesh::NODE,AmanziMesh::OWNED);
  unsigned int num_cells = mesh_maps.num_entities(AmanziMesh::CELL,AmanziMesh::OWNED);
	
  // create zone
  isize[0][0] = num_nodes;
  isize[1][0] = num_cells;
  isize[2][0] = 0;
  // goto file->base->end or use existing values
  cg_zone_write(file_idx, base_idx, "MeshZone", isize[0], Unstructured, &zone_idx);
	
  // get coords
  double *x = new double [num_nodes];
  double *y = new double [num_nodes];
  double *z = new double [num_nodes];
	
  for (int i=0; i<num_nodes; i++) {
    AmanziGeometry::Point xc(3);
    mesh_maps.node_get_coordinates(i,&xc);
    x[i] = xc[0];
    y[i] = xc[1];
    z[i] = xc[2];
  }
	
  // write out coords
  cg_coord_write(file_idx, base_idx, zone_idx, RealDouble, "CoordinateX", x, &coord_idx);
  cg_coord_write(file_idx, base_idx, zone_idx, RealDouble, "CoordinateY", y, &coord_idx);
  cg_coord_write(file_idx, base_idx, zone_idx, RealDouble, "CoordinateZ", z, &coord_idx);
	
  // get connectivity
  int *ielem = new int [num_cells*8];
  AmanziMesh::Entity_ID_List xh(8);
	
  nelem_start = 1; // first cell ID
  nelem_end = 0;
  for (int i=0; i<num_cells; i++) {
    mesh_maps.cell_get_nodes(i,&xh);
    for (int j=0; j<xh.size(); j++) {
      ielem[i*8+j] = xh[j]+1;
    }
    nelem_end++;
  }
  nbdyelem    = 0;		// unsorted boundary elements
	
  // write out connectivity
  cg_section_write(file_idx, base_idx, zone_idx, "Cells", HEXA_8, \
                   nelem_start, nelem_end, nbdyelem, ielem, &section_idx);
	
  // free memory
  delete x;
  delete y;
  delete z;
  delete ielem;
	
  // close file 
  cg_close(file_idx);
}
    
void open_data_file( std::string mesh_filename,  std::string soln_filename, 
                     unsigned int num_nodes, unsigned int num_cells)
{
  int isize[3][1];
	
  // open file
  //cg_set_file_type(CG_FILE_ADF);
  cg_open((const char*)soln_filename.c_str(), CG_MODE_WRITE, &file_idx);
	
  // create base - does it already exist?
  int icelldim = 3;  // cell dimension: 3-volume, 2-face, ...
  int iphysdim = 3;  // physical dimension: 1-1D, 2-2D, 3-3D
  cg_base_write(file_idx, "Base", icelldim, iphysdim, &base_idx);
  cg_simulation_type_write(file_idx, base_idx, TimeAccurate);
	
  // create zone through link nodes, need meshfilename, link path (/Base/zonename) - already exist?
  // create zone
  isize[0][0] = num_nodes;
  isize[1][0] = num_cells;
  isize[2][0] = 0;
  cg_zone_write(file_idx, base_idx, "MeshZone", isize[0], Unstructured, &zone_idx);
  //cg_goto(file_idx, base_idx, zone_idx,"end");
  /*
    cg_link_write("CoordinateX", (const char*)mesh_filename.c_str(), "/Base/MeshZone/CoordinateX");
    cg_link_write("CoordinateY", (const char*)mesh_filename.c_str(), "/Base/MeshZone/CoordinateY");
    cg_link_write("CoordinateZ", (const char*)mesh_filename.c_str(), "/Base/MeshZone/CoordinateZ");
    cg_link_write("Cells", (const char*)mesh_filename.c_str(), "/Base/MeshZone/Cells");
  */
  cg_link_write("GridCoordinates", (const char*)mesh_filename.c_str(), "/Base/MeshZone/GridCoordinates");
  cg_link_write("Cells", (const char*)mesh_filename.c_str(), "/Base/MeshZone/Cells");
	 
	
}
    
void open_data_file(const std::string soln_filename)
{
  // open file
  cg_open((const char*)soln_filename.c_str(), CG_MODE_MODIFY, &file_idx);
  base_idx = 1;
  zone_idx = 1;
	
}
    
void create_timestep(const double time, const int iter, AmanziMesh::Entity_kind kind)
{
  int size[2];
  std::stringstream solname;
		
  // make timestep name, add to list of names
  solname << "Solution" << iter;
  nameList.push_back(solname.str());
	
  // add current time to timeList
  timeList.push_back(time);
  num_soln++;
	
  // reconstruct formatted string of names to pass to cgns
  char solutionList[32*num_soln+1];
  sprintf(solutionList,"%-32s",nameList[0].c_str());
  for (int i=1; i<nameList.size(); i++) {
    sprintf(solutionList,"%s%-32s",solutionList,nameList[i].c_str());
  }
	
  // EIB: CGNS does not take advantage of HDF5's ability to append datasets; 
  //      therefore, we must track the time values and solution names, and overwrite
  //      these tables every time we add a solution.
  //      Also, CGNS is in C, so the list of solution names is taken as a 
  //      C char of size [32*numNames]
	
  // rewrite base iteration
  cg_biter_write(file_idx, base_idx,"TimeIterValues",num_soln);
  cg_goto(file_idx, base_idx,"BaseIterativeData_t",1,"end");
  cg_array_write("TimeValues",RealDouble,1,&num_soln,&timeList[0]);
	
  // rewrite zone iteration
  cg_ziter_write(file_idx, base_idx,zone_idx,"ZoneIterativeData");
  cg_goto(file_idx, base_idx,"Zone_t",zone_idx,"ZoneIterativeData_t",1,"end");
  size[0] = 32;
  size[1] = num_soln;
  cg_array_write("SolutionPointers",Character,2,size,solutionList);
	
  // create solution node under zone
  if (kind == AmanziMesh::CELL) {
    cg_sol_write(file_idx, base_idx,zone_idx, nameList.back().c_str(), CellCenter, &soln_idx);
  } else if (kind == AmanziMesh::NODE) {
    cg_sol_write(file_idx, base_idx,zone_idx, nameList.back().c_str(), Vertex, &soln_idx);
  } else {
    //throw an error
    printf("  E>> SHOULD BE THROWING ERROR HERE\n");
  }

	
}
    
void write_field_data(const Epetra_Vector &x, std::string varname)
{
  // write field data
  double *data;
  int err = x.ExtractView(&data);
  cg_field_write(file_idx, base_idx,zone_idx, soln_idx, RealDouble, (char *)varname.c_str(), data, &field_idx);
	
}
   
    
/*
  void write_probe_data(const Epetra_Vector &x, const Epetra_Vector &data)
  {
  }
*/
    
void close_data_file()
{
  cg_close(file_idx);
}
    
} // close namespace CGNS

} // close namespace Amanzi
