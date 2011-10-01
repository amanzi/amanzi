#ifndef __Transport_State_hpp__
#define __Transport_State_hpp__

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"


namespace Amanzi {
namespace AmanziTransport {

enum TransportCreateMode {
CopyPointers,  // copy Teuchos::RCP pointers 
ViewMemory,    // convert to overlap to non-overlap vectors  
CopyMemory     // copy non-overlap vector to overlap vectors 
};


/* The transport state is equivalent to the global state. */
class Transport_State {
 public:
  Transport_State() {};
  explicit Transport_State(State& S);
  Transport_State(Transport_State& S, TransportCreateMode mode = CopyPointers);
  ~Transport_State() {};

  // data management
  void copymemory_multivector(Epetra_MultiVector& source, Epetra_MultiVector& target);
  void copymemory_vector(Epetra_Vector& source, Epetra_Vector& target);

  // access methods for state variables
  Teuchos::RCP<Epetra_MultiVector> get_total_component_concentration() { return total_component_concentration; }

  Teuchos::RCP<Epetra_Vector> get_porosity() { return porosity; }
  Teuchos::RCP<Epetra_Vector> get_water_saturation() { return water_saturation; }
  Teuchos::RCP<Epetra_Vector> get_darcy_flux() { return darcy_flux; }
  Teuchos::RCP<Epetra_Vector> get_water_density() { return water_density; }
  Teuchos::RCP<AmanziMesh::Mesh> get_mesh_maps() { return mesh_maps; }

  Epetra_MultiVector& ref_total_component_concentration() { return *total_component_concentration; }

  Epetra_Vector& ref_porosity() { return *porosity; }
  Epetra_Vector& ref_water_saturation() { return *water_saturation; }
  Epetra_Vector& ref_darcy_flux() { return *darcy_flux; }
  Epetra_Vector& ref_water_density() { return *water_density; }
  
  // debug routines
  void analytic_total_component_concentration(double f(double*, double), double t = 0.0);
  void analytic_total_component_concentration(double tcc);
  void analytic_darcy_flux(double* u);
  void analytic_porosity(double phi = 0.2);
  void analytic_water_saturation(double ws = 1.0);
  void analytic_water_density(double wd = 1000.0);

  void error_total_component_concentration(double f(double*, double), double t, double* L1, double* L2);

 private:
  // state variables that are relevant to transport 
  Teuchos::RCP<Epetra_MultiVector> total_component_concentration;
  Teuchos::RCP<Epetra_Vector> water_saturation;
  Teuchos::RCP<Epetra_Vector> darcy_flux;
  Teuchos::RCP<Epetra_Vector> porosity;
  Teuchos::RCP<Epetra_Vector> water_density;

  Teuchos::RCP<AmanziMesh::Mesh> mesh_maps;
};

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif

