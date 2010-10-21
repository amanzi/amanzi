#include "State.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Mesh_maps_simple.hh"
extern "C" {
#include "gmvwrite.h"
}


State::State( int number_of_components_,
	      Teuchos::RCP<Mesh_maps_simple> mesh_maps_):
  number_of_components(number_of_components_),
  mesh_maps(mesh_maps_)
{
  // create the Eptera_Vector objects

  water_density =    Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  pressure =         Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  darcy_flux =       Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  porosity =         Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) );
  water_saturation = Teuchos::rcp( new Epetra_Vector( mesh_maps->cell_map(false) ) ); 
  total_component_concentration 
     = Teuchos::rcp( new Epetra_MultiVector( mesh_maps->cell_map(false), number_of_components ) );  

};


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

  double xc[3];
  for (int i=0; i<num_nodes; i++) {
    mesh_maps->node_to_coordinates(i,xc,xc+3);

    x[i] = xc[0];
    y[i] = xc[1];
    z[i] = xc[2];
  }
  gmvwrite_node_data(&num_nodes, x, y, z);

  delete x;
  delete y;
  delete z;
  
  unsigned int num_cells = mesh_maps->count_entities(Mesh_data::CELL,OWNED);
  
  gmvwrite_cell_header(&num_cells);
  
  int *xh = new int[8]; 
  for (int i=0; i<num_cells; i++) {
    mesh_maps->cell_to_nodes(i,xh,xh+8);
    for (int j=0; j<8; j++) xh[j]++;
    gmvwrite_cell_type("phex8",8,xh);
  }

  
  // write the side sets
  unsigned int num_side_sets = mesh_maps->num_sets(Mesh_data::FACE);
  if (num_side_sets > 0) {
    vector<unsigned int> ssids(num_side_sets);
    mesh_maps->get_set_ids(Mesh_data::FACE, ssids.begin(), ssids.end());
    

    vector<unsigned int> side_set;
    int total_num_surfaces=0;
    for (vector<unsigned int>::const_iterator it = ssids.begin(); it != ssids.end(); it++) {
      unsigned int set_size = mesh_maps->get_set_size(*it, Mesh_data::FACE, OWNED);
      total_num_surfaces += set_size;
    }
    if (total_num_surfaces > 0) {
      gmvwrite_surface_header(&total_num_surfaces);
      
      for (vector<unsigned int>::const_iterator it = ssids.begin(); it != ssids.end(); it++) {  
	unsigned int set_size = mesh_maps->get_set_size(*it, Mesh_data::FACE, OWNED);
	side_set.resize(set_size);
	
	mesh_maps->get_set(*it,Mesh_data::FACE,OWNED,side_set.begin(),side_set.end());
	
	unsigned int nodes[4];
	for (vector<unsigned int>::const_iterator sit = side_set.begin(); sit != side_set.end(); sit++) {
	  mesh_maps->face_to_nodes(*sit,nodes,nodes+4);
	  for (int j=0; j<4; j++) nodes[j]++;
	  gmvwrite_surface_data(4,nodes);
	}
      }
      
      
      // for now this only works with Mesh_maps_simple which has 
      // exactly six side sets
      gmvwrite_surfflag_header();
      gmvwrite_surfflag_name("SideSets",num_side_sets);
      gmvwrite_surfflag_subname("00");
      gmvwrite_surfflag_subname("01");
      gmvwrite_surfflag_subname("02");
      gmvwrite_surfflag_subname("03");
      gmvwrite_surfflag_subname("04");
      gmvwrite_surfflag_subname("05");
      
      int *sflagdata = new int[total_num_surfaces];
      unsigned int is=0;
      unsigned int tag=1;
      for (vector<unsigned int>::const_iterator it = ssids.begin(); it != ssids.end(); it++) {  
	unsigned int set_size = mesh_maps->get_set_size(*it, Mesh_data::FACE, OWNED);
	side_set.resize(set_size);
	
	for (int i=0; i<set_size; i++) {
	  sflagdata[is] = tag;
	  is++;
	}
	tag++;
      }
      gmvwrite_surfflag_data(sflagdata);

      gmvwrite_surfflag_endflag();  

    }
    
  }
  
  gmvwrite_closefile();
}
