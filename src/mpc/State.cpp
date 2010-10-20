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
    cout << xc[0] << " " << xc[1] << " " << xc[2] << endl;

    x[i] = xc[0];
    y[i] = xc[1];
    z[i] = xc[2];
  }

  unsigned int num_cells = mesh_maps->count_entities(Mesh_data::CELL,OWNED);

  gmvwrite_node_data(&num_nodes, x, y, z);
  
  gmvwrite_cell_header(&num_cells);
  
  int *xh = new int[8]; 
  for (int i=0; i<num_cells; i++) {
    mesh_maps->cell_to_nodes(i,xh,xh+8);
    for (int j=0; j<8; j++) xh[j]++;
    gmvwrite_cell_type("hex",8,xh);
  }

  gmvwrite_closefile();
}
