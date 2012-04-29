#ifndef __Transport_State_hpp__
#define __Transport_State_hpp__

#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"


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
  void copymemory_multivector(Epetra_MultiVector& source, Epetra_MultiVector& target, int target_is_parallel = 1);
  void copymemory_vector(Epetra_Vector& source, Epetra_Vector& target);
  void distribute_cell_vector(Epetra_Vector& v);
  void distribute_cell_multivector(Epetra_MultiVector& v);

  // access methods for state variables
  Teuchos::RCP<Epetra_MultiVector> total_component_concentration() { return total_component_concentration_; }

  Teuchos::RCP<Epetra_Vector> porosity() { return porosity_; }
  Teuchos::RCP<Epetra_Vector> water_saturation() { return water_saturation_; }
  Teuchos::RCP<Epetra_Vector> prev_water_saturation() { return prev_water_saturation_; }  
  Teuchos::RCP<Epetra_Vector> darcy_flux() { return darcy_flux_; }
  Teuchos::RCP<Epetra_Vector> water_density() { return water_density_; }
  Teuchos::RCP<AmanziMesh::Mesh> mesh() { return mesh_; }

  Epetra_MultiVector& ref_total_component_concentration() { return *total_component_concentration_; }

  Epetra_Vector& ref_porosity() { return *porosity_; }
  Epetra_Vector& ref_water_saturation() { return *water_saturation_; }
  Epetra_Vector& ref_prev_water_saturation() { return *prev_water_saturation_; }  
  Epetra_Vector& ref_darcy_flux() { return *darcy_flux_; }
  Epetra_Vector& ref_water_density() { return *water_density_; }

  // miscaleneous
  void interpolateCellVector(
      const Epetra_Vector& v0, const Epetra_Vector& v1, double dT_int, double dT, Epetra_Vector& v_int);
  double get_time() { return S_->get_time(); }
  
  // debug routines
  void analytic_total_component_concentration(double f_tcc(const AmanziGeometry::Point&, double), double t = 0.0);
  void analytic_total_component_concentration(double tcc);
  void analytic_darcy_flux(const AmanziGeometry::Point& u);
  void analytic_darcy_flux(AmanziGeometry::Point f_vel(const AmanziGeometry::Point&, double), double t = 0.0);
  void analytic_porosity(double phi = 0.2);
  void analytic_water_saturation(double ws = 1.0);
  void analytic_water_density(double wd = 1000.0);

  void error_total_component_concentration(double f_tcc(const AmanziGeometry::Point&, double), double t, double* L1, double* L2);

 private:
  State* S_;  
  Teuchos::RCP<Epetra_MultiVector> total_component_concentration_;
  Teuchos::RCP<Epetra_Vector> water_saturation_;
  Teuchos::RCP<Epetra_Vector> prev_water_saturation_;
  Teuchos::RCP<Epetra_Vector> darcy_flux_;
  Teuchos::RCP<Epetra_Vector> porosity_;
  Teuchos::RCP<Epetra_Vector> water_density_;

  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
};

}  // namespace AmanziTransport
}  // namespace Amanzi

#endif

