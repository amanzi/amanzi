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
  explicit Transport_State(State_Old& S);
  Transport_State(Transport_State& S, TransportCreateMode mode = CopyPointers);
  ~Transport_State() {};

  // data management
  void CopyMasterCell2GhostCell(Epetra_Vector& v);
  void CopyMasterCell2GhostCell(const Epetra_Vector& v, Epetra_Vector& vv);
  void CopyMasterFace2GhostFace(const Epetra_Vector& v, Epetra_Vector& vv);
  void CopyMasterMultiCell2GhostMultiCell(Epetra_MultiVector& v);
  void CopyMasterMultiCell2GhostMultiCell(const Epetra_MultiVector& v, 
                                          Epetra_MultiVector& vv, int parallel_comm = 1);

  // extension of Trilinos
  void MinValueMasterCells(Epetra_MultiVector& v, double* vmin);
  void MaxValueMasterCells(Epetra_MultiVector& v, double* vmax);

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
  void InterpolateCellVector(
      const Epetra_Vector& v0, const Epetra_Vector& v1, double dT_int, double dT, Epetra_Vector& v_int);
  double intermediate_time() { return (S_ == NULL) ? -1.0 : S_->intermediate_time(); }
  double initial_time() { return (S_ == NULL) ? -1.0 : S_->initial_time(); }
  double final_time() { return (S_ == NULL) ? -1.0 : S_->final_time(); }
  
  // debug routines
  void AnalyticTotalComponentConcentration(double f_tcc(const AmanziGeometry::Point&, double), double t = 0.0);
  void AnalyticTotalComponentConcentration(double tcc);
  void AnalyticDarcyFlux(const AmanziGeometry::Point& u);
  void AnalyticDarcyFlux(AmanziGeometry::Point f_vel(const AmanziGeometry::Point&, double), double t = 0.0);
  void AnalyticPorosity(double phi = 0.2);
  void AnalyticWaterSaturation(double ws = 1.0);
  void AnalyticWaterDensity(double wd = 1000.0);

  void error_total_component_concentration(double f_tcc(const AmanziGeometry::Point&, double), double t, double* L1, double* L2);

  // component name and number accessors
  int get_component_number(const std::string component_name) { return S_->get_component_number(component_name); }
  std::string get_component_name(const int component_number) { return S_->get_component_name(component_number); }

 private:
  State_Old* S_;  
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

