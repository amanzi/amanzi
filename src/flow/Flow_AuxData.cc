#include "Darcy_PK.hh"
#include "Richards_PK.hh"


namespace Amanzi {
namespace AmanziFlow {

void aux_compute_hydraulic_head(Epetra_Vector* hydraulic_head, double p_atm, Epetra_Vector* pressure, double rho,
                                AmanziGeometry::Point gravity, Epetra_Vector* centroids) {
  int dim = gravity.dim();

  hydraulic_head->PutScalar(-p_atm);
  hydraulic_head->Update(1.0,*pressure,1.0);
  
  double g = gravity[dim-1];
  
  hydraulic_head->Scale(1.0/(g*rho));
  hydraulic_head->Update(1.0, *centroids, 1.0);
}



void Darcy_PK::UpdateAuxilliaryData() {
  // update hydraulic head
  
  Epetra_Vector& pressure = FS->ref_pressure();
  Epetra_Vector& hydraulic_head = FS->ref_hydraulic_head();

  Epetra_Vector z_centroid(pressure);
  
  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);  
  int dim = FS->gravity()->dim();

  for (int c = 0; c != ncells; ++c) {
    const AmanziGeometry::Point& xp = mesh_->cell_centroid(c); 
    z_centroid[c] = xp[dim-1]; 
  }
  double rho = *FS->fluid_density();
  const AmanziGeometry::Point gravity = *FS->gravity();

  aux_compute_hydraulic_head(&hydraulic_head, atm_pressure, &pressure, rho, gravity,
                             &z_centroid);
}



void Richards_PK::UpdateAuxilliaryData() {
  // update hydraulic head
  std::cout << "UpdateAuxilliaryData\n";  
  Epetra_Vector& pressure = FS->ref_pressure();
  Epetra_Vector& hydraulic_head = FS->ref_hydraulic_head();
  
  Epetra_Vector z_centroid(pressure);

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);  
  int dim = FS->gravity()->dim();

  for (int c = 0; c != ncells; ++c) {
    const AmanziGeometry::Point& xp = mesh_->cell_centroid(c); 
    z_centroid[c] = xp[dim-1]; 
  }
  double rho = *FS->fluid_density();
  const AmanziGeometry::Point gravity = *FS->gravity();

  aux_compute_hydraulic_head(&hydraulic_head, atm_pressure, &pressure, rho, gravity,
                             &z_centroid);
}




}
}
