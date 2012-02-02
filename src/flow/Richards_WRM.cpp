/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (nnc@lanl.gov), 
         Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "Richards_PK.hpp"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Defines relative permeability ONLY for cells.                                               
****************************************************************** */
void Richards_PK::calculateRelativePermeabilityCell(const Epetra_Vector& p)
{
  for (int mb=0; mb<WRM.size(); mb++) {
    std::string region = WRM[mb]->region();

    //AmanziMesh::Set_ID_List block;
    std::vector<AmanziMesh::Set_ID> block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    std::vector<unsigned int>::iterator i;
    for (i=block.begin(); i!=block.end(); i++) {
      double pc = atm_pressure - p[*i];
      (*Krel_cells)[*i] = WRM[mb]->k_relative(pc);
    }
  }
}


/* ******************************************************************
* Wrapper for various ways to define relative permeability of faces.
****************************************************************** */
void Richards_PK::calculateRelativePermeabilityFace(const Epetra_Vector& p)
{
  calculateRelativePermeabilityCell(p);  // populates cell-based permeabilities
  FS->copyMasterCell2GhostCell(*Krel_cells);

  if (Krel_method == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) {  // Define K and Krel_faces
    calculateRelativePermeabilityUpwindGravity(p);
  } 
  else if (Krel_method == FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX) {
    Epetra_Vector& flux = FS_next->ref_darcy_mass_flux();

    matrix->createMFDstiffnessMatrices(mfd3d_method, K, *Krel_faces);
    matrix->deriveDarcyFlux(*solution, *face_importer_, flux);
    addGravityFluxes_DarcyFlux(K, *Krel_faces, flux);

    calculateRelativePermeabilityUpwindFlux(p, flux);
  } 
  else if (Krel_method == FLOW_RELATIVE_PERM_ARITHMETIC_MEAN) {
    calculateRelativePermeabilityArithmeticMean(p);
  }
}


/* ******************************************************************
* Defines upwinded relative permeabilities for faces using gravity. 
****************************************************************** */
void Richards_PK::calculateRelativePermeabilityUpwindGravity(const Epetra_Vector& p)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c=0; c<ncells_owned; c++) {
    mesh_->cell_get_face_dirs(c, &dirs);
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    AmanziGeometry::Point Kgravity = K[c] * gravity;

    for (int n=0; n<nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f); 
      if ((normal * Kgravity) * dirs[n] >= 0.0) (*Krel_faces)[f] = (*Krel_cells)[c];
      else if (bc_markers[f] != FLOW_BC_FACE_NULL) (*Krel_faces)[f] = (*Krel_cells)[c];
    }
  }
}


/* ******************************************************************
* Defines upwinded relative permeabilities for faces using a given flux. 
****************************************************************** */
void Richards_PK::calculateRelativePermeabilityUpwindFlux(const Epetra_Vector& p, 
                                                          const Epetra_Vector& flux)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c=0; c<ncells_owned; c++) {
    mesh_->cell_get_face_dirs(c, &dirs);
    mesh_->cell_get_faces(c, &faces);
    int nfaces = faces.size();

    for (int n=0; n<nfaces; n++) {
      int f = faces[n];
      if (flux[n] * dirs[n] >= 0.0) (*Krel_faces)[f] = (*Krel_cells)[c];
      else if (bc_markers[f] != FLOW_BC_FACE_NULL) (*Krel_faces)[f] = (*Krel_cells)[c];
    }
  }
}


/* ******************************************************************
* Defines relative permeabilities for faces via arithmetic averaging. 
****************************************************************** */
void Richards_PK::calculateRelativePermeabilityArithmeticMean(const Epetra_Vector& p)
{
  AmanziMesh::Entity_ID_List cells;

  Krel_faces->PutScalar(0.0);
  for (int f=0; f<nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    for (int n=0; n<ncells; n++) (*Krel_faces)[f] += (*Krel_cells)[cells[n]];
    (*Krel_faces)[f] /= ncells;
  }
}


/* ******************************************************************
* Use analytical formula for derivative dS/dP.                                               
****************************************************************** */
void Richards_PK::derivedSdP(const Epetra_Vector& p, Epetra_Vector& ds)
{
  for (int mb=0; mb<WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    int ncells = mesh_->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);

    std::vector<unsigned int> block(ncells);
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);
      
    std::vector<unsigned int>::iterator i;
    for (i=block.begin(); i!=block.end(); i++) {
      double pc = atm_pressure - p[*i];
      ds[*i] = WRM[mb]->d_saturation(pc);
    }
  }
}


/* ******************************************************************
* Convertion p -> s.                                               
****************************************************************** */
void Richards_PK::deriveSaturationFromPressure(const Epetra_Vector& p, Epetra_Vector& s)
{
  for (int mb=0; mb<WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    int ncells = mesh_->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);

    std::vector<unsigned int> block(ncells);
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);
      
    std::vector<unsigned int>::iterator i;
    for (i=block.begin(); i!=block.end(); i++) {
      double pc = atm_pressure - p[*i];
      s[*i] = WRM[mb]->saturation(pc);
    }
  }
}


/* ******************************************************************
* Convertion s -> p.                                               
****************************************************************** */
void Richards_PK::derivePressureFromSaturation(double s, Epetra_Vector& p)
{
  for (int mb=0; mb<WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    int ncells = mesh_->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);

    std::vector<unsigned int> block(ncells);
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);
      
    std::vector<unsigned int>::iterator i;
    for (i=block.begin(); i!=block.end(); i++) {
      double pc = WRM[mb]->capillaryPressure(s);
      p[*i] = atm_pressure - pc;
    }
  } 
 
  for (int mb=0; mb<WRM.size(); mb++) {
    std::string region = WRM[mb]->region();
    int ncells = mesh_->get_set_size(region, AmanziMesh::CELL, AmanziMesh::OWNED);

    std::vector<unsigned int> block(ncells);
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);
      
    std::vector<unsigned int>::iterator i;
    for (i=block.begin(); i!=block.end(); i++) {
      double pc = WRM[mb]->capillaryPressure(s);
      p[*i] = atm_pressure - pc;
    }
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi



