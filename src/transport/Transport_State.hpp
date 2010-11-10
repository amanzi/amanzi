#ifndef __Transport_State_hpp__
#define __Transport_State_hpp__

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

#include "Mesh_maps_base.hh"


/* 
  The transport state is a sub-state of the global state.
  At the moment they are equivalent. 
*/


using namespace Teuchos;


class Transport_State {

public:
  Transport_State ( State S );
  /* a null constructor is useful for unit tests */
  Transport_State () {};

  /* major member functions */
  void copy_constant_state( Transport_State & TS );
  void create_internal_state( Transport_State & TS );

  ~Transport_State () {};

  /* access methods for state variables */
  RCP<Epetra_MultiVector>   get_total_component_concentration() const { return total_component_concentration; }
  RCP<Epetra_MultiVector>   get_total_component_concentration()       { return total_component_concentration; }

  RCP<const Epetra_Vector>  get_porosity() const { return porosity; }
  RCP<Epetra_Vector>        get_porosity()       { return porosity; }

  RCP<const Epetra_Vector>  get_water_saturation() const { return water_saturation; }
  RCP<Epetra_Vector>        get_water_saturation()       { return water_saturation; }

  RCP<const Epetra_Vector>  get_darcy_flux() const { return darcy_flux; }
  RCP<Epetra_Vector>        get_darcy_flux()       { return darcy_flux; }
  
  RCP<const Epetra_Vector>  get_water_density() const { return water_density; }
  RCP<Epetra_Vector>        get_water_density()       { return water_density; }
  
  RCP<Mesh_maps_base> get_mesh_maps() const { return mesh_maps; }

  /* debug routines */
  void analytic_total_component_concentration();
  void analytic_darcy_flux( double* u );
  void analytic_porosity( double phi = 0.2 );
  void analytic_water_saturation( double ws = 1.0 );
  void analytic_water_density( double wd = 1000.0 );


private:
  /* state variables that are relevant to transport */
  RCP<Epetra_MultiVector>   total_component_concentration;
  RCP<Epetra_Vector>        water_saturation;
  RCP<Epetra_Vector>        darcy_flux;
  RCP<Epetra_Vector>        porosity;
  RCP<Epetra_Vector>        water_density;

  /* mesh infranstructure */
  RCP<Mesh_maps_base> mesh_maps;
};


#endif

