#ifndef __Status_hpp__
#define __Status_hpp__

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "DataLayout.hpp"
#include "MeshWrapper.hpp"


typedef enum { COMPLETE, UPDATING } status_type;

class State {

public:

  State( Teuchos::RCP<DataLayout>, int, Teuchos::RCP<MeshWrapper> );
  ~State() {};

  // access methods

  Teuchos::RCP<const Epetra_Vector> get_pressure () const { return pressure; }; 
  Teuchos::RCP<const Epetra_Vector> get_darcy_flux () const { return darcy_flux; };
  Teuchos::RCP<const Epetra_Vector> get_porosity () const { return porosity; };
  Teuchos::RCP<const Epetra_Vector> get_water_saturation () const { return water_saturation; };
  Teuchos::RCP<const Epetra_Vector> get_water_density () const { return water_density; };
  Teuchos::RCP<Epetra_MultiVector> get_total_component_concentration () 
  { return total_component_concentration; };
  
  Teuchos::RCP<const MeshWrapper> get_mesh_wrapper() const { return mesh_wrapper; };

  const double get_time () const { return time; };

  // modify methods

  void set_time ( double new_time );

  // status methods
  
  const status_type get_status () const { return status; };
  void set_status ( status_type new_status ) { status = new_status; }
      
private:
  Teuchos::RCP<Epetra_Vector> water_density;  
  Teuchos::RCP<Epetra_Vector> pressure;
  Teuchos::RCP<Epetra_Vector> darcy_flux;
  Teuchos::RCP<Epetra_Vector> porosity;
  Teuchos::RCP<Epetra_MultiVector> total_component_concentration; 
  Teuchos::RCP<Epetra_Vector> water_saturation;
  
  int number_of_components;

  double time;
  status_type status;

  // data layout
  Teuchos::RCP<DataLayout> data_layout;
  
  // mesh
  Teuchos::RCP<MeshWrapper> mesh_wrapper;
}; 


#endif
