#ifndef __State_hpp__
#define __State_hpp__

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "../simple_mesh/Mesh_maps_base.hh"


typedef enum { COMPLETE, UPDATING } status_type;


class State {

public:

  State( int, Teuchos::RCP<Mesh_maps_base> );

  State( Teuchos::ParameterList &parameter_list, 
	 Teuchos::RCP<Mesh_maps_base> );

  ~State() {};

  // access methods

  Teuchos::RCP<const Epetra_Vector> get_pressure () const { return pressure; }; 

  Teuchos::RCP<const Epetra_Vector> get_darcy_flux () const { return darcy_flux; };
  Teuchos::RCP<Epetra_Vector>       get_darcy_flux ()       { return darcy_flux; };

  Teuchos::RCP<const Epetra_Vector> get_porosity () const { return porosity; };
  Teuchos::RCP<Epetra_Vector>       get_porosity ()       { return porosity; };

  Teuchos::RCP<const Epetra_Vector> get_water_saturation () const { return water_saturation; };
  Teuchos::RCP<Epetra_Vector>       get_water_saturation ()       { return water_saturation; };

  Teuchos::RCP<const Epetra_Vector> get_water_density () const { return water_density; };
  Teuchos::RCP<Epetra_Vector>       get_water_density ()       { return water_density; };

  Teuchos::RCP<Epetra_MultiVector> get_total_component_concentration () 
  { return total_component_concentration; };
  
  const Teuchos::RCP<Mesh_maps_base> get_mesh_maps() const { return mesh_maps; };

  const double get_time () const { return time; };

  const int get_number_of_components() const { return number_of_components; };

  // modify methods
  void set_time ( double new_time );
  void advance_time(double dT);
  void update_total_component_concentration(Teuchos::RCP<Epetra_MultiVector>);

  // status methods
  
  const status_type get_status () const { return status; };
  void set_status ( status_type new_status ) { status = new_status; }


  void write_gmv ( std::string filename );
      
private:
  void create_storage();
  

  // state vectors
  Teuchos::RCP<Epetra_Vector> water_density;  
  Teuchos::RCP<Epetra_Vector> pressure;
  Teuchos::RCP<Epetra_Vector> darcy_flux;
  Teuchos::RCP<Epetra_Vector> porosity;
  Teuchos::RCP<Epetra_MultiVector> total_component_concentration; 
  Teuchos::RCP<Epetra_Vector> water_saturation;
  
  int number_of_components;

  double time;
  status_type status;

  // mesh
  const Teuchos::RCP<Mesh_maps_base> mesh_maps;
}; 


#endif
