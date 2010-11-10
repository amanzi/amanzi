#ifndef __Flow_State_hpp__
#define __Flow_State_hpp__

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "Mesh_maps_base.hh"
#include "State.hpp"

class Flow_State {

public:
    
  Flow_State(Teuchos::RCP<State> S):
    mesh_maps_(S->get_mesh_maps())
    // INITIALIZE THE OTHER DATA COMPONENTS
  { };

  ~Flow_State () {};

  // access methods
  const Teuchos::RCP<Mesh_maps_base>& mesh() const { return mesh_maps_;};
  
  double fluid_density () const { return *fluid_density_; }
  
  double fluid_viscosity () const { return *fluid_viscosity_; }
  
  const double* gravity() const { return *gravity_; }
  
  const std::vector<double>& permeability() const { return *permeability_; }

private:
    
  // object doesn't own anything -- all smart pointers to the real thing.
    
  const Teuchos::RCP<const double> fluid_density_;
  const Teuchos::RCP<const double> fluid_viscosity_;
  const Teuchos::RCP<const double[3]> gravity_;
  const Teuchos::RCP<const std::vector<double> > permeability_;
  //const Teuchos::RCP<const std::vector<Epetra_SerialSymDenseMatrix>> permeability_;
  const Teuchos::RCP<Mesh_maps_base> mesh_maps_;
};



#endif
