#include "State.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Mesh_maps_base.hh"
extern "C" {
#include "gmvwrite.h"
}


State::State( int number_of_components_,
	      Teuchos::RCP<Mesh_maps_base> mesh_maps_):
  number_of_components(number_of_components_),
  mesh_maps(mesh_maps_)
{
  // create the Eptera_Vector objects

  create_storage();

};

State::State( Teuchos::ParameterList &parameter_list,
	      Teuchos::RCP<Mesh_maps_base> mesh_maps_):
  mesh_maps(mesh_maps_)
{

  // get the number of component concentrations from the 
  // parameter list
  number_of_components = parameter_list.get<int>("Number of component concentrations");
  
  // create the Eptera_Vector objects

  create_storage();

};


void State::create_storage ()
{
  // create the Eptera_Vector objects

  water_density =    Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  pressure =         Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  darcy_flux =       Teuchos::rcp( new Epetra_Vector( mesh_maps->face_map(false) ) );
  porosity =         Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  water_saturation = Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) ); 
  total_component_concentration 
    = Teuchos::rcp( new Epetra_MultiVector( mesh_maps->cell_map(false), number_of_components ) );  

}



void State::set_time ( double new_time ) {

  if ( status == UPDATING ) {
    
    time = new_time;

  } else {
    
    // throw an error

  }

}


void State::write_gmv ( std::string filename )
{
  using namespace std;

  gmvwrite_openfile_ir_ascii( (char*) filename.c_str(), 4, 8 );
  
  // first write the node data
  unsigned int num_nodes = mesh_maps->count_entities(Mesh_data::NODE,OWNED);

  double *x = new double [num_nodes];
  double *y = new double [num_nodes];
  double *z = new double [num_nodes];

  std::vector<double> xc(3);
  for (unsigned int i=0; i < num_nodes; i++) {
    mesh_maps->node_to_coordinates(i,xc.begin(),xc.end());

    x[i] = xc[0];
    y[i] = xc[1];
    z[i] = xc[2];
  }
  gmvwrite_node_data(&num_nodes, x, y, z);

  delete x;
  delete y;
  delete z;
  

  // write the cell data
  unsigned int num_cells = mesh_maps->count_entities(Mesh_data::CELL,OWNED);
  
  gmvwrite_cell_header(&num_cells);
  
  std::vector<unsigned int> xh(8); 
  unsigned int xh_[8];
  for (int i=0; i<num_cells; i++) {
    mesh_maps->cell_to_nodes(i,xh.begin(),xh.end());
    for (int j=0; j<8; j++) xh_[j] = xh[j]+1;

    char cell_type [] = "phex8";
    
    
    gmvwrite_cell_type(cell_type,8,xh_);
  }
  

  
  // write the side sets, we only support 100 side sets
  // with ids less than 100.

  // how many side sets are there?
  unsigned int num_side_sets = mesh_maps->num_sets(Mesh_data::FACE);
  if (num_side_sets > 0) {
    vector<unsigned int> ssids(num_side_sets);
    mesh_maps->get_set_ids(Mesh_data::FACE, ssids.begin(), ssids.end());
    
    vector<unsigned int> side_set;
    int total_num_surfaces=0;

    // figure out how many total faces are in side sets
    for (vector<unsigned int>::const_iterator it = ssids.begin(); it != ssids.end(); it++) {
      unsigned int set_size = mesh_maps->get_set_size(*it, Mesh_data::FACE, OWNED);
      total_num_surfaces += set_size;
    }

    if (total_num_surfaces > 0) {
      // now start writing side set information
      gmvwrite_surface_header(&total_num_surfaces);
      
      // loop over the side sets and write all faces belonging to side sets
      for (vector<unsigned int>::const_iterator it = ssids.begin(); it != ssids.end(); it++) {  
	unsigned int set_size = mesh_maps->get_set_size(*it, Mesh_data::FACE, OWNED);
	side_set.resize(set_size);
	
	mesh_maps->get_set(*it,Mesh_data::FACE,OWNED,side_set.begin(),side_set.end());
	
	std::vector<unsigned int> nodes(4);
	unsigned int nodes_[4];
	for (vector<unsigned int>::const_iterator sit = side_set.begin(); sit != side_set.end(); sit++) {
	  mesh_maps->face_to_nodes(*sit,nodes.begin(),nodes.end());
	  for (int j=0; j<4; j++) nodes_[j] = nodes[j]+1;
	  gmvwrite_surface_data(4,nodes_);
	}
      }
      
      
      // now we generate an array of flags, one flag stands for one side set
      gmvwrite_surfflag_header();
      
      char ssname [] = "SideSets";
      gmvwrite_surfflag_name(ssname,num_side_sets);
      
      char sssubname [] = "SideSet_000";

      for (vector<unsigned int>::const_iterator it=ssids.begin(); 
	   it!=ssids.end(); it++) {
	sssubname[8] = '0' + (*it)/100;
	sssubname[9] = '0' + ((*it)/10)%10;
	sssubname[10] = '0' + (*it)%10;

	gmvwrite_surfflag_subname(sssubname);
      }
      
      
      // this is the array of flags
      int *sflagdata = new int[total_num_surfaces];
      unsigned int is=0;
      unsigned int tag=1;
      
      // loop over the side sets
      for (vector<unsigned int>::const_iterator it = ssids.begin(); it != ssids.end(); it++) {  
	unsigned int set_size = mesh_maps->get_set_size(*it, Mesh_data::FACE, OWNED);
	side_set.resize(set_size);
	
	// loop over the faces in the current side set and
	// store the flag for the current side set 
 	for (int i=0; i<set_size; i++) {
	  sflagdata[is] = tag;
	  is++;
	}
	tag++;
      }
      // write the flag array
      gmvwrite_surfflag_data(sflagdata);
      
      delete sflagdata;
      gmvwrite_surfflag_endflag();
    }
    
  }
  

  // write element blocks

  unsigned int num_element_blocks = mesh_maps->num_sets(Mesh_data::CELL);
  if (num_element_blocks > 0) {
    vector<unsigned int> ebids(num_element_blocks);
    mesh_maps->get_set_ids(Mesh_data::CELL, ebids.begin(), ebids.end());
    
    vector<unsigned int> element_block;
    int total_num_cells=0;

    // figure out how many total cells are in element blocks
    for (vector<unsigned int>::const_iterator it = ebids.begin(); it != ebids.end(); it++) {
      unsigned int set_size = mesh_maps->get_set_size(*it, Mesh_data::CELL, OWNED);
      total_num_cells += set_size;
    }

    if (total_num_cells > 0) {
      // now start writing side set information
      gmvwrite_flag_header();

      char flagname [] = "ElementBlocks";
      gmvwrite_flag_name(flagname,num_element_blocks,0);


      // now we generate an array of flags, one flag stands for one side set
      char ebsubname [] = "ElementBlock_000";

      for (vector<unsigned int>::const_iterator it=ebids.begin(); 
	   it!=ebids.end(); it++) {
	ebsubname[13] = '0' + (*it)/100;
	ebsubname[14] = '0' + ((*it)/10)%10;
	ebsubname[15] = '0' + (*it)%10;

	gmvwrite_flag_subname(ebsubname);
      }
      
      // loop over the element blocks and write all cells belonging to element blocks

      // this is the array of flags
      int *ebflagdata = new int[total_num_cells];
      unsigned int is=0;
      unsigned int tag=1;
      
      // loop over the selement blocks sets
      for (vector<unsigned int>::const_iterator it = ebids.begin(); it != ebids.end(); it++) {  
	unsigned int set_size = mesh_maps->get_set_size(*it, Mesh_data::CELL, OWNED);
	
	// loop over the cells in the current element block set and
	// store the flag for the current side set 
 	for (int i=0; i<set_size; i++) {
	  ebflagdata[is] = tag;
	  is++;
	}
	tag++;
      }
      // write the flag array
      gmvwrite_flag_data(0,ebflagdata);

      delete ebflagdata;

      gmvwrite_flag_endflag();
    }
    
  }  
 

  // done
  gmvwrite_closefile();
}
