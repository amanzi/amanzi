/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Neil Carlson (nnc@lanl.gov), 
         Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <string>
#include <vector>

#include "Richards_PK.hh"


namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Wrapper for various ways to define dKdP on faces.
****************************************************************** */
void Richards_PK::CalculateDerivativePermeabilityFace(const Epetra_Vector& p)
{
  DerivedKdP(p, *dKdP_cells);  // populates cell-based permeabilities
  FS->CopyMasterCell2GhostCell(*dKdP_cells);

  if (Krel_method == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) {
    CalculateDerivativePermeabilityUpwindGravity(p);
  } else if (Krel_method == FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX) {
    Epetra_Vector& flux = FS_aux->ref_darcy_flux();
    CalculateDerivativeRelativePermeabilityUpwindFlux(p, flux);
  }
}


/* ******************************************************************
* Defines upwind value of dKdP on faces using gravity. 
****************************************************************** */
void Richards_PK::CalculateDerivativePermeabilityUpwindGravity(const Epetra_Vector& p)
{
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm->WRM();  

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  dKdP_faces->PutScalar(0.0);

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    AmanziGeometry::Point Kg = K[c] * gravity_;
    Kg /= norm(Kg);

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      double cos_angle = (normal * Kg) * dirs[n] / mesh_->face_area(f);
      
      if (bc_model[f] != FLOW_BC_FACE_NULL){
        if (bc_model[f] == FLOW_BC_FACE_PRESSURE && 
            cos_angle < -FLOW_RELATIVE_PERM_TOLERANCE) {
          int mb = (rel_perm->map_c2mb())[c];
          double pc = atm_pressure - bc_values[f][0];
          (*dKdP_faces)[f] = WRM[mb]->dKdPc(pc);
        } else {
          (*dKdP_faces)[f] = (*dKdP_cells)[c];
        }
      } else {
        if (cos_angle > FLOW_RELATIVE_PERM_TOLERANCE) {
          (*dKdP_faces)[f] = (*dKdP_cells)[c];  // The upwind face.
        } else if (fabs(cos_angle) <= FLOW_RELATIVE_PERM_TOLERANCE) { 
          (*dKdP_faces)[f] += (*dKdP_cells)[c] / 2;  // Almost vertical face.
        }
      }
    }
  }
}

/* ******************************************************************
* Defines upwind derivative of relative permeability on mesh faces 
* using a given flux.
* WARNING: This is a part of the experimental solver. 
****************************************************************** */
void Richards_PK::CalculateDerivativeRelativePermeabilityUpwindFlux(
    const Epetra_Vector& p, const Epetra_Vector& flux)
{
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm->WRM();  

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  dKdP_faces->PutScalar(0.0);

  double max_flux;
  flux.MaxValue(&max_flux);
  double tol = FLOW_RELATIVE_PERM_TOLERANCE * max_flux;

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];

      if (bc_model[f] != FLOW_BC_FACE_NULL) {  // The boundary face.
        if ((bc_model[f] == FLOW_BC_FACE_PRESSURE) &&
            (flux[f] * dirs[n] < -tol)) {
          int mb = (rel_perm->map_c2mb())[c];
          double pc = atm_pressure - bc_values[f][0];
          (*dKdP_faces)[f] = WRM[mb]->dKdPc(pc);
        } else {
          (*dKdP_faces)[f] = (*dKdP_cells)[c];
        }
      } else {
        if (flux[f] * dirs[n] > tol) {
          (*dKdP_faces)[f] = (*dKdP_cells)[c];  // The upwind face.
        } else if (fabs(flux[f]) <= tol) { 
          (*dKdP_faces)[f] += (*dKdP_cells)[c] / 2; // Zero flux face.
        }
      }
    }
  }
}


/* ******************************************************************
* Use analytical formula for derivative dS/dP.                                               
****************************************************************** */
void Richards_PK::DerivedSdP(const Epetra_Vector& p, Epetra_Vector& ds)
{
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm->WRM();  

  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double pc = atm_pressure - p[*i];
      ds[*i] = -WRM[mb]->dSdPc(pc);  // Negative sign indicates that dSdP = -dSdPc.
    }
  }
}


/* ******************************************************************
* Use analytical formula for derivative dK/dP.                                               
****************************************************************** */
void Richards_PK::DerivedKdP(const Epetra_Vector& p, Epetra_Vector& dk)
{
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm->WRM();  

  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double pc = atm_pressure - p[*i];
      dk[*i] = -WRM[mb]->dKdPc(pc);  // Negative sign indicates that dKdP = -dKdPc.
    }
  }
}


/* ******************************************************************
* Convertion p -> s.                                               
****************************************************************** */
void Richards_PK::DeriveSaturationFromPressure(const Epetra_Vector& p, Epetra_Vector& s)
{
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm->WRM();  

  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double pc = atm_pressure - p[*i];
      s[*i] = WRM[mb]->saturation(pc);
    }
  }
}


/* ******************************************************************
* Convertion s -> p.                                               
****************************************************************** */
void Richards_PK::DerivePressureFromSaturation(const Epetra_Vector& s, Epetra_Vector& p)
{
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm->WRM();  

  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double pc = WRM[mb]->capillaryPressure(s[*i]);
      p[*i] = atm_pressure - pc;
    }
  }
}


/* ******************************************************************
* Clip pressure using pressure threshold.
****************************************************************** */
void Richards_PK::ClipHydrostaticPressure(const double pmin, Epetra_Vector& p)
{
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm->WRM();  

  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    // double pc = atm_pressure - pmin;
    // double s0 = WRM[mb]->saturation(pc);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      if (p[*i] < pmin) p[*i] = pmin;
    }
  }
}


/* ******************************************************************
* Clip pressure using constant saturation.
****************************************************************** */
void Richards_PK::ClipHydrostaticPressure(const double pmin, const double s0, Epetra_Vector& p)
{
  std::vector<Teuchos::RCP<WaterRetentionModel> >& WRM = rel_perm->WRM();  

  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      if (p[*i] < pmin) {
        double pc = WRM[mb]->capillaryPressure(s0);
        p[*i] = atm_pressure - pc;
      }
    }
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi



