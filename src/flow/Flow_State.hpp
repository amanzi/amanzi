#ifndef __Flow_State_hpp__
#define __Flow_State_hpp__

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "Mesh.hh"
#include "State.hpp"

namespace Amanzi {

class Flow_State {
 public:
  Flow_State(Teuchos::RCP<State> S) :
      mesh_maps_(S->get_mesh_maps()),
      gravity_(S->get_gravity()),
      fluid_density_(S->get_density()),
      fluid_viscosity_(S->get_viscosity()),
      vertical_permeability_(S->get_vertical_permeability()),
      horizontal_permeability_(S->get_horizontal_permeability()),
      pressure_(S->get_pressure()),
      porosity_(S->get_porosity()),
      water_saturation_(S->get_water_saturation()),
      prev_water_saturation_(S->get_prev_water_saturation()) 
  {};
  
  ~Flow_State() {};

  // access methods
  const Teuchos::RCP<AmanziMesh::Mesh>& get_mesh_maps() const { return mesh_maps_ ;}
  double get_fluid_density() const { return *fluid_density_; }
  double get_fluid_viscosity() const { return *fluid_viscosity_; }
  const double* get_gravity() const { return *gravity_; }
  const Epetra_Vector& get_vertical_permeability() const { return *vertical_permeability_; }
  const Epetra_Vector& get_horizontal_permeability() const { return *horizontal_permeability_; }
  const Epetra_Vector& get_porosity() const { return *porosity_; }
  Epetra_Vector& get_water_saturation() { return *water_saturation_; }
  Epetra_Vector& get_prev_water_saturation() { return *prev_water_saturation_; }
  Epetra_Vector& get_pressure() { return *pressure_; }

 private:
  // object doesn't own anything -- all smart pointers to the real thing.
  const Teuchos::RCP<double> fluid_density_;
  const Teuchos::RCP<double> fluid_viscosity_;
  const Teuchos::RCP<double*> gravity_;
  const Teuchos::RCP<Epetra_Vector> vertical_permeability_;
  const Teuchos::RCP<Epetra_Vector> horizontal_permeability_;
  const Teuchos::RCP<AmanziMesh::Mesh> mesh_maps_;
  const Teuchos::RCP<Epetra_Vector> pressure_;  // current cell pressure solution
  const Teuchos::RCP<Epetra_Vector> porosity_;
  const Teuchos::RCP<Epetra_Vector> water_saturation_;
  const Teuchos::RCP<Epetra_Vector> prev_water_saturation_;  
};

} // close namespace Amanzi

#endif
