#include "cgns_mesh.hh"
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_IntVector.h"
//#include "Entity_kind.hh"
//#include "Element_category.hh"
//#include <vector>



namespace CGNS_PAR {

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
  int rank;
  std::vector<double> timeList;
  std::vector<int> cycleList;
  std::vector<string> nameList;
  Teuchos::RCP<Epetra_Map> all_to_one_node_map;
  Teuchos::RCP<Epetra_Map> all_to_one_cell_map;  
  Teuchos::RCP<Epetra_Export> all_to_one_node_export;
  Teuchos::RCP<Epetra_Export> all_to_one_cell_export;
  

  
  
  void create_mesh_file(Mesh_maps_base &mesh_maps, std::string filename)
  {
    using namespace CGNS_PAR;
    
    
    int isize[3][1];
    int nelem_start, nelem_end, nbdyelem;
    
    rank = mesh_maps.get_comm()->MyPID();
    
    // open file
    if (rank == 0 ) cg_open ((const char*)filename.c_str(), CG_MODE_WRITE, &file_idx);
    
    // create base
    int icelldim = 3;  // cell dimension: 3-volume, 2-face, ...
    int iphysdim = 3;  // physical dimension: 1-1D, 2-2D, 3-3D
    
    if (rank==0) cg_base_write(file_idx, "Base", icelldim, iphysdim, &base_idx);
    
    // get num_nodes, num_cells
    // - gather sum to PE0 - MB
    unsigned int num_nodes = mesh_maps.count_entities(Mesh_data::NODE,OWNED);
    unsigned int num_cells = mesh_maps.count_entities(Mesh_data::CELL,OWNED);
    
    int nums[2];
    int dummy[2];
    dummy[0] = num_nodes;
    dummy[1] = num_cells;
    
    
    mesh_maps.get_comm()->SumAll(dummy,nums,2);

    if (rank==0) {
      // create zone
      isize[0][0] = nums[0];
      isize[1][0] = nums[1];
      isize[2][0] = 0;
      // goto file->base->end or use existing values
      cg_zone_write(file_idx, base_idx, "MeshZone", isize[0], Unstructured, &zone_idx);
    }
    
    // get coords
    // - instead use Epetra_Vector for coordinates and communicate to 
    // - PE0 with export
    double *x = new double [num_nodes];
    double *y = new double [num_nodes];
    double *z = new double [num_nodes];
    
    std::vector<double> xc(3);
    for (int i=0; i<num_nodes; i++) {
      mesh_maps.node_to_coordinates(i,xc.begin(),xc.end());
      x[i] = xc[0];
      y[i] = xc[1];
      z[i] = xc[2];
    }
    // make Epetra versions of the coordinate vectors
    Epetra_Vector ep_x(View, mesh_maps.node_map(false), x);
    Epetra_Vector ep_y(View, mesh_maps.node_map(false), y);
    Epetra_Vector ep_z(View, mesh_maps.node_map(false), z);
    
    // make the all to one map
    if (rank == 0) {
      int *gids = new int[nums[1]];
      for (int i=0; i<nums[1]; i++) gids[i] = i;
      all_to_one_cell_map = Teuchos::rcp(new Epetra_Map(nums[1],nums[1],gids,0, *mesh_maps.get_comm() ));
      delete [] gids;
      
      gids = new int[nums[0]];
      for (int i=0; i<nums[0]; i++) gids[i] = i;    
      all_to_one_node_map = Teuchos::rcp(new Epetra_Map(nums[0],nums[0],gids,0, *mesh_maps.get_comm() ));
      delete [] gids;
      
    } else {
      int *gids;
      all_to_one_cell_map = Teuchos::rcp(new Epetra_Map(nums[1],0,gids,0,*mesh_maps.get_comm() ) );
      all_to_one_node_map = Teuchos::rcp(new Epetra_Map(nums[0],0,gids,0,*mesh_maps.get_comm() ) );
    }
    
    //all_to_one_cell_map->Print(cout);
    //(mesh_maps.cell_map(false)).Print(cout); 

    // make the all to one exporters 
    all_to_one_cell_export = Teuchos::rcp(new Epetra_Export(mesh_maps.cell_map(false), *all_to_one_cell_map) );
    all_to_one_node_export = Teuchos::rcp(new Epetra_Export(mesh_maps.node_map(false), *all_to_one_node_map) );
    
    // transfer x coordinates to PE0
    Epetra_Vector PE0(*all_to_one_node_map);
    
    double **view = new double*;
    
    // gather x coordinates 
    PE0.Export(ep_x, *all_to_one_node_export, Insert);
    // - pass the view of the Epetra_Vector PE0
    PE0.ExtractView(view);
    // write out x coordinates
    if (rank == 0) cg_coord_write(file_idx, base_idx, zone_idx, RealDouble, "CoordinateX", *view, &coord_idx);
    
    // gather y coordinates 
    PE0.Export(ep_y, *all_to_one_node_export, Insert);
    // - pass the view of the Epetra_Vector PE0
    PE0.ExtractView(view);
    // write out y coordinates
    if (rank == 0) cg_coord_write(file_idx, base_idx, zone_idx, RealDouble, "CoordinateY", *view, &coord_idx);
    
    // gather z coordinates 
    PE0.Export(ep_z, *all_to_one_node_export, Insert);
    // - pass the view of the Epetra_Vector PE0
    PE0.ExtractView(view);
    // write out z coordinates
    if (rank == 0) cg_coord_write(file_idx, base_idx, zone_idx, RealDouble, "CoordinateZ", *view, &coord_idx);
    
    // get connectivity
    int *ielem = new int [num_cells*8];
    std::vector<unsigned int> xh(8);
    
    nelem_start = 1; // first cell ID
    nelem_end = 0;
    for (int i=0; i<num_cells; i++) {
      mesh_maps.cell_to_nodes(i,xh.begin(),xh.end());
      for (int j=0; j<8; j++) {
	ielem[i*8+j] = mesh_maps.node_map(true).GID(xh[j]) +1;
      }
      nelem_end++;
    }
    nbdyelem    = 0;		// unsorted boundary elements
    
    int *gids = new int [num_cells*8];
    int LIDStart = mesh_maps.cell_map(false).MinLID();
    for (int i=0; i<num_cells; i++) {
      int start =  8*mesh_maps.cell_map(false).GID(LIDStart+i);
      for (int j=0; j<8; j++) gids[8*i+j] = start + j;
    }
    
    Epetra_Map cell_to_node_map(nums[1]*8,num_cells*8,gids,0,*mesh_maps.get_comm());
    delete [] gids;
    //cell_to_node_map.Print(cout);

    // make an epetra vector of the cell to node data
    Epetra_IntVector ivec(View, cell_to_node_map, ielem);
    //ivec.Print(cout);
    
    // create an all to one map for the cell to node data
    int num_loc;
    if (rank == 0) {
      gids = new int[nums[1]*8];
      for (int i=0; i<8*nums[1]; i++) gids[i] = i;
      num_loc = 8*nums[1];
    } else {
      num_loc = 0;
    }   
    Epetra_Map all_to_one_cell_to_nodes_map(8*nums[1],num_loc,gids,0,*mesh_maps.get_comm());
    if (!gids) delete [] gids;
  
    //all_to_one_cell_to_nodes_map.Print(cout); 

    //cout << rank << " here"<< endl;
  
    // create an all to one exporter for the cell to node data
    Epetra_Export all_to_one_cell_to_nodes_export(cell_to_node_map, all_to_one_cell_to_nodes_map);
    //all_to_one_cell_to_nodes_export.Print(cout);

    //cout << rank << " here 1"<< endl;

    // transfer the cell to node data to PE0
    Epetra_IntVector iPE0(all_to_one_cell_to_nodes_map);
    iPE0.Export(ivec,all_to_one_cell_to_nodes_export,Insert);
    //iPE0.Print(cout);

    // write out connectivity
    // - write on PE0 using a view of an Epetra_IntVector
    
    int **iview = new int*;
    iPE0.ExtractView(iview);
    
    if (rank==0) cg_section_write(file_idx, base_idx, zone_idx, "Cells", HEXA_8, 
    				  nelem_start, nums[1], nbdyelem, *iview, &section_idx);
    
    // free memory
    delete x;
    delete y;
    delete z;
    delete ielem;
    
    // close file 
    if (rank == 0) cg_close(file_idx);
  }
  
  void open_data_file( std::string mesh_filename,  std::string soln_filename, 
		       unsigned int num_nodes, unsigned int num_cells)
  {
    using namespace CGNS_PAR;
    int isize[3][1];

    // open file
    //cg_set_file_type(CG_FILE_ADF);
    if (rank == 0) { 
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
      
  }
    
  void open_data_file(const std::string soln_filename)
  {
    using namespace CGNS_PAR;

    if (rank == 0) {
      // open file
      cg_open((const char*)soln_filename.c_str(), CG_MODE_MODIFY, &file_idx);
      base_idx = 1;
      zone_idx = 1;
      
    }
  }
    
  void create_timestep(const double time, const int iter, Mesh_data::Entity_kind kind)
    {
      using namespace CGNS_PAR;

      if (rank == 0) {
	
	int size[2];
	std::stringstream solname;
	
	// make timestep name, add to list of names
	solname << "Solution" << iter;
	nameList.push_back(solname.str());
	
	// add current time to timeList
	timeList.push_back(time);
	num_soln++;

	// add the current cycle to the list of cycles
	cycleList.push_back(iter);
		
	
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
	cg_array_write("Cycles",Integer,1,&num_soln,&cycleList[0]);

	// rewrite zone iteration
	cg_ziter_write(file_idx, base_idx,zone_idx,"ZoneIterativeData");
	cg_goto(file_idx, base_idx,"Zone_t",zone_idx,"ZoneIterativeData_t",1,"end");
	size[0] = 32;
	size[1] = num_soln;
	cg_array_write("SolutionPointers",Character,2,size,solutionList);
	
	
	
	// create solution node under zone
	if (kind == Mesh_data::CELL) {
	  cg_sol_write(file_idx, base_idx,zone_idx, nameList.back().c_str(), CellCenter, &soln_idx);
	} else if (kind == Mesh_data::NODE) {
	  cg_sol_write(file_idx, base_idx,zone_idx, nameList.back().c_str(), Vertex, &soln_idx);
	} else {
	  //throw an error
	  printf("  E>> SHOULD BE THROWING ERROR HERE\n");
	}
	
      }
    }
    
  void write_field_data(const Epetra_Vector &x, std::string varname)
  {
    using namespace CGNS_PAR;
    
    Epetra_Vector PE0(*all_to_one_cell_map);    
    PE0.Export(x, *all_to_one_cell_export,Insert);
    
    if (rank == 0) {
      
      // write field data
      double *data;
      // - export the vector to PE0 and write
      int err = PE0.ExtractView(&data);
      cg_field_write(file_idx, base_idx,zone_idx, soln_idx, RealDouble, (char *)varname.c_str(), data, &field_idx);
    }
  }
  
    
  /*
    void write_probe_data(const Epetra_Vector &x, const Epetra_Vector &data)
    {
    }
  */
  
  void close_data_file()
  {
    using namespace CGNS_PAR;
    
    if (rank == 0) cg_close(file_idx);
  }
  
}
