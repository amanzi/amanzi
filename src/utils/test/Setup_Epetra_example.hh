#ifndef _SETUP_EPETRA_EXAMPLE_H_
#define _SETUP_EPETRA_EXAMPLE_H_

//#include "test/Setup_tests.hh"
#include <Teuchos_RCP.hpp>
#include "Epetra_Vector.h"

using namespace std;

typedef Teuchos::RCP<Epetra_Vector> EVec_p;

struct Epetra_mesh_setup : public Map_setup
{
  struct Map_setup stk_map_setup;
  STK_mesh::Mesh_maps_stk stk_map;
  EVec_p node_quantity;
  EVec_p cell_quantity;
  EVec_p fake_pressure;

  Epetra_mesh_setup() :
    stk_map_setup(),
    stk_map(stk_map_setup.mesh_map)
  {
    // Setup node quantity
    int node_index_list[] = {5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    double node_values[] = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120};
    node_quantity = Teuchos::rcp( new Epetra_Vector(stk_map.node_map(false)));
    node_quantity->ReplaceGlobalValues(12, node_values, node_index_list);

    // Setup cell quantity
    int cell_index_list[] = {1, 2, 3, 4};
    double cell_values[] = {10, 20, 30, 40};
    cell_quantity = Teuchos::rcp( new Epetra_Vector(stk_map.cell_map(false))); 
    cell_quantity->ReplaceGlobalValues(4, cell_values, cell_index_list);

    // Setup second cell quantity -- called fake pressure
    double fake_values[] = {9, 8, 7, 6};
    fake_pressure = Teuchos::rcp( new Epetra_Vector(stk_map.cell_map(false))); 
    fake_pressure->ReplaceGlobalValues(4, fake_values, cell_index_list);
 }

};


#endif /* _SETUP_EPETRA_EXAMPLE_H_ */
