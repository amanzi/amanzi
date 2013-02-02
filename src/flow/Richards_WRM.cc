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
* Defines relative permeability ONLY for cells.                                               
****************************************************************** */
void Richards_PK::CalculateRelativePermeabilityCell(const Epetra_Vector& p)
{
  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();

    std::vector<AmanziMesh::Entity_ID> block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double pc = atm_pressure - p[*i];
      (*Krel_cells)[*i] = WRM[mb]->k_relative(pc);
    }
  }
}


/* ******************************************************************
* Wrapper for various ways to define relative permeability of faces.
****************************************************************** */
void Richards_PK::CalculateRelativePermeabilityFace(const Epetra_Vector& p)
{
  CalculateRelativePermeabilityCell(p);  // populates cell-based permeabilities
  FS->CopyMasterCell2GhostCell(*Krel_cells);

  if (Krel_method == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) {  // Define K and Krel_faces
    CalculateRelativePermeabilityUpwindGravity(p);

  } else if (Krel_method == FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX) {
    Epetra_Vector& flux = FS_aux->ref_darcy_flux();
    // const Epetra_Map& fmap = mesh_->face_map(true);
    // Epetra_Vector flux_wghost(fmap);
    // FS->CopyMasterFace2GhostFace(flux, flux_wghost);
    CalculateRelativePermeabilityUpwindFlux(p, flux);

  } else if (Krel_method == FLOW_RELATIVE_PERM_ARITHMETIC_MEAN) {
    CalculateRelativePermeabilityArithmeticMean(p);

  } else if (Krel_method == FLOW_RELATIVE_PERM_EXPERIMENTAL) {
    CalculateRelativePermeabilityUpwindGravity(p);
    AverageRelativePermeability();
  }
}


/* ******************************************************************
* Defines upwinded relative permeabilities for faces using gravity. 
****************************************************************** */
void Richards_PK::CalculateRelativePermeabilityUpwindGravity(const Epetra_Vector& p)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Krel_faces->PutScalar(0.0);

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      double cos_angle = (normal * Kgravity_unit[c]) * dirs[n] / mesh_->face_area(f);

      if (bc_model[f] != FLOW_BC_FACE_NULL) {  // The boundary face.
        if (bc_model[f] == FLOW_BC_FACE_PRESSURE && cos_angle < -FLOW_RELATIVE_PERM_TOLERANCE) {
          double pc = atm_pressure - bc_values[f][0];
          (*Krel_faces)[f] = WRM[(*map_c2mb)[c]]->k_relative(pc);
        } else {
          (*Krel_faces)[f] = (*Krel_cells)[c];
        }
      } else {
        if (cos_angle > FLOW_RELATIVE_PERM_TOLERANCE) {
          (*Krel_faces)[f] = (*Krel_cells)[c]; // The upwind face.
        } else if (fabs(cos_angle) <= FLOW_RELATIVE_PERM_TOLERANCE) { 
          (*Krel_faces)[f] += (*Krel_cells)[c] / 2;  //Almost vertical face.
        } 
      }
    }
  }
}


/* ******************************************************************
* Defines upwinded relative permeabilities for faces using a given flux.
* WARNING: This is the experimental code. 
****************************************************************** */
void Richards_PK::CalculateRelativePermeabilityUpwindFlux(const Epetra_Vector& p,
                                                          const Epetra_Vector& flux)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Krel_faces->PutScalar(0.0);

  double max_flux, min_flux;
  flux.MaxValue(&max_flux);
  flux.MinValue(&min_flux);
  double tol = FLOW_RELATIVE_PERM_TOLERANCE * std::max<double>(fabs(max_flux), fabs(min_flux));

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];

      if (bc_model[f] != FLOW_BC_FACE_NULL) {  // The boundary face.
        if (bc_model[f] == FLOW_BC_FACE_PRESSURE && flux[f] * dirs[n] < -tol) {
          double pc = atm_pressure - bc_values[f][0];
          (*Krel_faces)[f] = WRM[(*map_c2mb)[c]]->k_relative(pc);
        } else {
          (*Krel_faces)[f] = (*Krel_cells)[c];
        }
      } else {
        if (flux[f] * dirs[n] > tol) {
          (*Krel_faces)[f] = (*Krel_cells)[c];  // The upwind face.
        } else if (fabs(flux[f]) <= tol) { 
          (*Krel_faces)[f] += (*Krel_cells)[c] / 2;  // Almost vertical face.
        }
      }
    }
  }
}


/* ******************************************************************
* Defines relative permeabilities for faces via arithmetic averaging. 
****************************************************************** */
void Richards_PK::CalculateRelativePermeabilityArithmeticMean(const Epetra_Vector& p)
{
  AmanziMesh::Entity_ID_List cells;

  Krel_faces->PutScalar(0.0);
  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    for (int n = 0; n < ncells; n++) (*Krel_faces)[f] += (*Krel_cells)[cells[n]];
    (*Krel_faces)[f] /= ncells;
  }
}


/* ******************************************************************
* Defines upwinded relative permeabilities for faces using gravity. 
****************************************************************** */
void Richards_PK::AverageRelativePermeability()
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    double factor = 0.0;
    for (int n = 0; n < nfaces; n++) factor += (*Krel_faces)[faces[n]];
    (*Krel_cells)[c] = factor / nfaces;
  } 
}


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
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  dKdP_faces->PutScalar(0.0);

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      double cos_angle = (normal * Kgravity_unit[c]) * dirs[n] / mesh_->face_area(f);
      
      if (bc_model[f] != FLOW_BC_FACE_NULL){
        if (bc_model[f] == FLOW_BC_FACE_PRESSURE && 
            cos_angle < -FLOW_RELATIVE_PERM_TOLERANCE) {
          int mb = (*map_c2mb)[c];
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
          int mb = (*map_c2mb)[c];
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


/* ******************************************************************
* Calculates tensor-point product K * g and normalizes the result.
* To minimize parallel communications, the resultin vector Kg_unit 
* is distributed across mesh.
****************************************************************** */
void Richards_PK::CalculateKVectorUnit(const AmanziGeometry::Point& g, 
                                       std::vector<AmanziGeometry::Point>& Kg_unit)
{
  const Epetra_Map& cmap = mesh_->cell_map(true);
  Epetra_MultiVector Kg_copy(cmap, dim);  // temporary vector

  for (int c = 0; c < ncells_owned; c++) {
    AmanziGeometry::Point Kg = K[c] * g;
    double Kg_norm = norm(Kg);
 
    for (int i = 0; i < dim; i++) Kg_copy[i][c] = Kg[i] / Kg_norm;
  }

#ifdef HAVE_MPI
  FS->CopyMasterMultiCell2GhostMultiCell(Kg_copy);
#endif

  Kg_unit.clear();
  for (int c = 0; c < ncells_wghost; c++) {
    AmanziGeometry::Point Kg(dim); 
    for (int i = 0; i < dim; i++) Kg[i] = Kg_copy[i][c];
    Kg_unit.push_back(Kg);
  }
} 


/* ******************************************************************
* Auxiliary map from cells to WRM models.                                               
****************************************************************** */
void Richards_PK::PopulateMapC2MB()
{
  map_c2mb->PutScalar(-1);

  for (int mb = 0; mb < WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) (*map_c2mb)[*i] = mb;
  }

  for (int c = 0; c < ncells_owned; c++) {
    if ((*map_c2mb)[c] < 0) {
      Errors::Message msg;
      msg << "Flow PK: water retention models do not cover the whole domain.";
      Exceptions::amanzi_throw(msg);  
    }
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi



