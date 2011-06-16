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
      permeability_(S->get_permeability()),
      pressure_(S->get_pressure()),
      porosity_(S->get_porosity())
    {}

  ~Flow_State () {};

  // access methods
  const Teuchos::RCP<AmanziMesh::Mesh>& mesh() const { return mesh_maps_;};

  double fluid_density () const { return *fluid_density_; }

  double fluid_viscosity () const { return *fluid_viscosity_; }

  const double* gravity() const { return *gravity_; }

  const Epetra_Vector& permeability() const { return *permeability_; }

  const Epetra_Vector& porosity() const { return *porosity_; }

private:

  // object doesn't own anything -- all smart pointers to the real thing.

  const Teuchos::RCP<double> fluid_density_;
  const Teuchos::RCP<double> fluid_viscosity_;
  const Teuchos::RCP<double*> gravity_;
  const Teuchos::RCP<Epetra_Vector> permeability_;
  const Teuchos::RCP<AmanziMesh::Mesh> mesh_maps_;
  const Teuchos::RCP<Epetra_Vector> pressure_;  // current cell pressure solution
  const Teuchos::RCP<Epetra_Vector> porosity_;
};

} // close namespace Amanzi

#endif
