/*
This is the flow component of the Amanzi code. 

Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "errors.hh"

#include "WaterRetentionModel.hh"
#include "WRM_vanGenuchten.hh"
#include "WRM_BrooksCorey.hh"
#include "WRM_fake.hh"

#include "RelativePermeability.hh"
#include "Flow_constants.hh"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Initialize internal data.
****************************************************************** */
void RelativePermeability::Init(double p0, const Teuchos::RCP<Flow_State> FS)
{
  atm_pressure = p0;
  FS_ = FS;

  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  Krel_cells_ = Teuchos::rcp(new Epetra_Vector(cmap_wghost));
  Krel_faces_ = Teuchos::rcp(new Epetra_Vector(fmap_wghost));
  dKdP_cells_ = Teuchos::rcp(new Epetra_Vector(cmap_wghost));
  dKdP_faces_ = Teuchos::rcp(new Epetra_Vector(fmap_wghost));

  upwind_cell = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));
  downwind_cell = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));
  face_flag = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));

  method_ = FLOW_RELATIVE_PERM_NONE;
  SetFullySaturated();

  map_c2mb_ = Teuchos::rcp(new Epetra_Vector(cmap_wghost));
}


/* ******************************************************************
* A wrapper for updating relative permeabilities.
****************************************************************** */
void RelativePermeability::Compute(const Epetra_Vector& p,
                                    const std::vector<int>& bc_model,
                                    const std::vector<bc_tuple>& bc_values)
{
  if (method_ == FLOW_RELATIVE_PERM_UPWIND_GRAVITY ||
      method_ == FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX ||
      method_ == FLOW_RELATIVE_PERM_ARITHMETIC_MEAN) {
    ComputeOnFaces(p, bc_model, bc_values);
    Krel_cells_->PutScalar(1.0);
    if (experimental_solver_ == FLOW_SOLVER_NEWTON || 
        experimental_solver_ == FLOW_SOLVER_PICARD_NEWTON) {
      ComputeDerivativeOnFaces(p, bc_model, bc_values);
    }
  } else if (method_ == FLOW_RELATIVE_PERM_AMANZI) {
    ComputeOnFaces(p, bc_model, bc_values);
  } else {
    ComputeInCells(p);
    Krel_faces_->PutScalar(1.0);
  }
}

 
/* ******************************************************************
* Defines relative permeability ONLY for cells.                                               
****************************************************************** */
void RelativePermeability::ComputeInCells(const Epetra_Vector& p)
{
  for (int mb = 0; mb < WRM_.size(); mb++) {
    std::string region = WRM_[mb]->region();

    std::vector<AmanziMesh::Entity_ID> block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double pc = atm_pressure - p[*i];
      (*Krel_cells_)[*i] = WRM_[mb]->k_relative(pc);
    }
  }
}


/* ******************************************************************
* Wrapper for various ways to define relative permeability of faces.
****************************************************************** */
void RelativePermeability::ComputeOnFaces(
    const Epetra_Vector& p,
    const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values)
{
  ComputeInCells(p);  // populates cell-based permeabilities
  FS_->CopyMasterCell2GhostCell(*Krel_cells_);

  if (method_ == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) {  // Define K and Krel_faces
    FaceUpwindGravity_(p, bc_model, bc_values);

  } else if (method_ == FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX) {
    Epetra_Vector& flux = FS_->ref_darcy_flux();
    FaceUpwindFlux_(p, flux, bc_model, bc_values);

  } else if (method_ == FLOW_RELATIVE_PERM_ARITHMETIC_MEAN) {
    FaceArithmeticMean_(p);

  } else if (method_ == FLOW_RELATIVE_PERM_AMANZI) {
    FaceUpwindGravityInSoil_(p, bc_model, bc_values);
  }
}


/* ******************************************************************
* Defines upwinded relative permeabilities for faces using gravity. 
****************************************************************** */
void RelativePermeability::FaceUpwindGravity_(
    const Epetra_Vector& p,
    const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Krel_faces_->PutScalar(0.0);

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      double cos_angle = (normal * Kgravity_unit_[c]) * dirs[n] / mesh_->face_area(f);

      if (bc_model[f] != FLOW_BC_FACE_NULL) {  // The boundary face.
        if (bc_model[f] == FLOW_BC_FACE_PRESSURE && cos_angle < -FLOW_RELATIVE_PERM_TOLERANCE) {
          double pc = atm_pressure - bc_values[f][0];
          (*Krel_faces_)[f] = WRM_[(*map_c2mb_)[c]]->k_relative(pc);
        } else {
          (*Krel_faces_)[f] = (*Krel_cells_)[c];
        }
      } else {
        if (cos_angle > FLOW_RELATIVE_PERM_TOLERANCE) {
          (*Krel_faces_)[f] = (*Krel_cells_)[c]; // The upwind face.
        } else if (fabs(cos_angle) <= FLOW_RELATIVE_PERM_TOLERANCE) { 
          (*Krel_faces_)[f] += (*Krel_cells_)[c] / 2;  //Almost vertical face.
        } 
      }
    }
  }
}


/* ******************************************************************
* Defines upwinded relative permeabilities for faces using gravity. 
****************************************************************** */
void RelativePermeability::FaceUpwindGravityInSoil_(
    const Epetra_Vector& p,
    const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Krel_amanzi_.clear();

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    std::vector<double> krel(nfaces);

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      int c1 = (*upwind_cell)[f];
      int c2 = (*downwind_cell)[f];

      int flag = (*face_flag)[f];
      if (flag == FLOW_PERMFLAG_AVERAGE) {
        krel[n] = ((*Krel_cells_)[c1] + (*Krel_cells_)[c2]) / 2;
        (*Krel_faces_)[f] = krel[n];
      } else if (flag == FLOW_PERMFLAG_INTERFACE) {
        krel[n] = (*Krel_cells_)[c];
        (*Krel_faces_)[f] = (*Krel_cells_)[c1];
      } else if (flag == FLOW_PERMFLAG_UPWIND) {
        if (bc_model[f] == FLOW_BC_FACE_PRESSURE && c1 < 0) {
          double pc = atm_pressure - bc_values[f][0];
          krel[n] = WRM_[(*map_c2mb_)[c]]->k_relative(pc);
        } else if (c1 >= 0) {
          krel[n] = (*Krel_cells_)[c1];
        } else {
          krel[n] = (*Krel_cells_)[c];
        }
        (*Krel_faces_)[f] = krel[n];
      } 
    }
    Krel_amanzi_.push_back(krel);
  }
}


/* ******************************************************************
* Defines upwinded relative permeabilities for faces using a given flux.
* WARNING: This is the experimental code. 
****************************************************************** */
void RelativePermeability::FaceUpwindFlux_(
    const Epetra_Vector& p, const Epetra_Vector& flux,
    const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Krel_faces_->PutScalar(0.0);

  double max_flux, min_flux;
  flux.MaxValue(&max_flux);
  flux.MinValue(&min_flux);
  double tol = FLOW_RELATIVE_PERM_TOLERANCE * std::max(fabs(max_flux), fabs(min_flux));

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];

      if (bc_model[f] != FLOW_BC_FACE_NULL) {  // The boundary face.
        if (bc_model[f] == FLOW_BC_FACE_PRESSURE && flux[f] * dirs[n] < -tol) {
          double pc = atm_pressure - bc_values[f][0];
          (*Krel_faces_)[f] = WRM_[(*map_c2mb_)[c]]->k_relative(pc);
        } else {
          (*Krel_faces_)[f] = (*Krel_cells_)[c];
        }
      } else {
        if (flux[f] * dirs[n] > tol) {
          (*Krel_faces_)[f] = (*Krel_cells_)[c];  // The upwind face.
        } else if (fabs(flux[f]) <= tol) { 
          (*Krel_faces_)[f] += (*Krel_cells_)[c] / 2;  // Almost vertical face.
        }
      }
    }
  }
}


/* ******************************************************************
* Defines relative permeabilities for faces via arithmetic averaging. 
****************************************************************** */
void RelativePermeability::FaceArithmeticMean_(const Epetra_Vector& p)
{
  AmanziMesh::Entity_ID_List cells;

  Krel_faces_->PutScalar(0.0);

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    for (int n = 0; n < ncells; n++) (*Krel_faces_)[f] += (*Krel_cells_)[cells[n]];
    (*Krel_faces_)[f] /= ncells;
  }
}


/* ******************************************************************
* Wrapper for various ways to define dKdP on faces.
****************************************************************** */
void RelativePermeability::ComputeDerivativeOnFaces(
    const Epetra_Vector& p,
    const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values)
{
  DerivedKdP(p, *dKdP_cells_);  // populates cell-based permeabilities
  FS_->CopyMasterCell2GhostCell(*dKdP_cells_);

  if (method_ == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) {
    DerivativeFaceUpwindGravity_(p, bc_model, bc_values);
  } else if (method_ == FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX) {
    Epetra_Vector& flux = FS_->ref_darcy_flux();
    DerivativeFaceUpwindFlux_(p, flux, bc_model, bc_values);
  }
}


/* ******************************************************************
* Defines upwind value of dKdP on faces using gravity. 
****************************************************************** */
void RelativePermeability::DerivativeFaceUpwindGravity_(
   const Epetra_Vector& p,
   const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  dKdP_faces_->PutScalar(0.0);

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      double cos_angle = (normal * Kgravity_unit_[c]) * dirs[n] / mesh_->face_area(f);
      
      if (bc_model[f] != FLOW_BC_FACE_NULL){
        if (bc_model[f] == FLOW_BC_FACE_PRESSURE && cos_angle < -FLOW_RELATIVE_PERM_TOLERANCE) {
          int mb = (*map_c2mb_)[c];
          double pc = atm_pressure - bc_values[f][0];
          (*dKdP_faces_)[f] = WRM_[mb]->dKdPc(pc);
        } else {
          (*dKdP_faces_)[f] = (*dKdP_cells_)[c];
        }
      } else {
        if (cos_angle > FLOW_RELATIVE_PERM_TOLERANCE) {
          (*dKdP_faces_)[f] = (*dKdP_cells_)[c];  // The upwind face.
        } else if (fabs(cos_angle) <= FLOW_RELATIVE_PERM_TOLERANCE) { 
          (*dKdP_faces_)[f] += (*dKdP_cells_)[c] / 2;  // Almost vertical face.
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
void RelativePermeability::DerivativeFaceUpwindFlux_(
    const Epetra_Vector& p, const Epetra_Vector& flux,
    const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  dKdP_faces_->PutScalar(0.0);

  double max_flux;
  flux.MaxValue(&max_flux);
  double tol = FLOW_RELATIVE_PERM_TOLERANCE * max_flux;

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];

      if (bc_model[f] != FLOW_BC_FACE_NULL) {  // The boundary face.
        if (bc_model[f] == FLOW_BC_FACE_PRESSURE && flux[f] * dirs[n] < -tol) {
          int mb = (*map_c2mb_)[c];
          double pc = atm_pressure - bc_values[f][0];
          (*dKdP_faces_)[f] = WRM_[mb]->dKdPc(pc);
        } else {
          (*dKdP_faces_)[f] = (*dKdP_cells_)[c];
        }
      } else {
        if (flux[f] * dirs[n] > tol) {
          (*dKdP_faces_)[f] = (*dKdP_cells_)[c];  // The upwind face.
        } else if (fabs(flux[f]) <= tol) { 
          (*dKdP_faces_)[f] += (*dKdP_cells_)[c] / 2; // Zero flux face.
        }
      }
    }
  }
}


/* ******************************************************************
* Use analytical formula for derivative dS/dP.                                               
****************************************************************** */
void RelativePermeability::DerivedSdP(const Epetra_Vector& p, Epetra_Vector& ds)
{
  for (int mb = 0; mb < WRM_.size(); mb++) {
    std::string region = WRM_[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double pc = atm_pressure - p[*i];
      ds[*i] = -WRM_[mb]->dSdPc(pc);  // Negative sign indicates that dSdP = -dSdPc.
    }
  }
}


/* ******************************************************************
* Use analytical formula for derivative dK/dP.                                               
****************************************************************** */
void RelativePermeability::DerivedKdP(const Epetra_Vector& p, Epetra_Vector& dk)
{
  for (int mb = 0; mb < WRM_.size(); mb++) {
    std::string region = WRM_[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double pc = atm_pressure - p[*i];
      dk[*i] = -WRM_[mb]->dKdPc(pc);  // Negative sign indicates that dKdP = -dKdPc.
    }
  }
}


/* ******************************************************************
* Defines upwinded relative permeabilities for faces using gravity.
* The unwind is restricted to each soil. 
****************************************************************** */
void RelativePermeability::FaceUpwindGravityInit_()
{
  for (int f = 0; f < nfaces_wghost; f++) {
    (*upwind_cell)[f] = -1;
    (*downwind_cell)[f] = -1;
    (*face_flag)[f] = FLOW_PERMFLAG_NONE;
  }

  // populate internal faces
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      double cos_angle = (normal * Kgravity_unit_[c]) * dirs[n] / mesh_->face_area(f);

      if (cos_angle < -FLOW_RELATIVE_PERM_TOLERANCE) {
        (*face_flag)[f] = FLOW_PERMFLAG_UPWIND;
        (*downwind_cell)[f] = c;  
      } else if (cos_angle > FLOW_RELATIVE_PERM_TOLERANCE) {
        (*face_flag)[f] = FLOW_PERMFLAG_UPWIND;
        (*upwind_cell)[f] = c;
      } else { 
        (*face_flag)[f] = FLOW_PERMFLAG_AVERAGE;  // Almost vertical face.
        if ((*upwind_cell)[f] < 0) 
          (*upwind_cell)[f] = c;  
        else
          (*downwind_cell)[f] = c; 
      }
    }
  }

  // update boundary faces
  for (int f = 0; f < nfaces_wghost; f++) {
    if ((*upwind_cell)[f] < 0 && (*face_flag)[f] == FLOW_PERMFLAG_AVERAGE)
       (*upwind_cell)[f] = (*downwind_cell)[f]; 
    if ((*downwind_cell)[f] < 0 && (*face_flag)[f] == FLOW_PERMFLAG_AVERAGE)
       (*downwind_cell)[f] = (*upwind_cell)[f]; 
  }

  // update internal interface faces
  for (int f = 0; f < nfaces_wghost; f++) {
    int c1 = (*upwind_cell)[f];
    int c2 = (*downwind_cell)[f];

    if (c1 >= 0 && c2 >= 0) {
      int mb1 = (*map_c2mb_)[c1];
      int mb2 = (*map_c2mb_)[c2];
      if (mb1 != mb2) (*face_flag)[f] = FLOW_PERMFLAG_INTERFACE;
    }
  }
}


/* *******************************************************************
* Identify flux direction based on the mutual orientation of the 
* face normal and gravity vector.  
* This is the obsolete routine.                      
******************************************************************* */
void RelativePermeability::FaceUpwindGravityInit_(const AmanziGeometry::Point& g)
{
  for (int f = 0; f < nfaces_owned; f++) {
    (*upwind_cell)[f] = -1;
    (*downwind_cell)[f] = -1;
  }

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> fdirs;

  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &fdirs);

    for (int i = 0; i < faces.size(); i++) {
      int f = faces[i];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      if ((normal * g) * fdirs[i] >= 0)
        (*upwind_cell)[f] = c;
      else
        (*downwind_cell)[f] = c;
    }
  }
}


/* ******************************************************************
* Calculates tensor-point product K * g and normalizes the result.
* To minimize parallel communications, the resultin vector Kg_unit 
* is distributed across mesh.
****************************************************************** */
void RelativePermeability::CalculateKVectorUnit(const std::vector<WhetStone::Tensor>& K,
                                                 const AmanziGeometry::Point& g)
{
  const Epetra_Map& cmap = mesh_->cell_map(true);
  int dim = g.dim();
  Epetra_MultiVector Kg_copy(cmap, dim);  // temporary vector

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells_owned; c++) {
    AmanziGeometry::Point Kg = K[c] * g;
    double Kg_norm = norm(Kg);
 
    for (int i = 0; i < dim; i++) Kg_copy[i][c] = Kg[i] / Kg_norm;
  }

#ifdef HAVE_MPI
  FS_->CopyMasterMultiCell2GhostMultiCell(Kg_copy);
#endif

  int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  Kgravity_unit_.clear();

  for (int c = 0; c < ncells_wghost; c++) {
    AmanziGeometry::Point Kg(dim); 
    for (int i = 0; i < dim; i++) Kg[i] = Kg_copy[i][c];
    Kgravity_unit_.push_back(Kg);
  }

  FaceUpwindGravityInit_();
} 


/* ******************************************************************
* Auxiliary map from cells to WRM models.                                               
****************************************************************** */
void RelativePermeability::PopulateMapC2MB()
{
  map_c2mb_->PutScalar(-1);

  for (int mb = 0; mb < WRM_.size(); mb++) {
    std::string region = WRM_[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) (*map_c2mb_)[*i] = mb;
  }
  
#ifdef HAVE_MPI
  FS_->CopyMasterCell2GhostCell(*map_c2mb_);
#endif

  // internal check
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells_owned; c++) {
    if ((*map_c2mb_)[c] < 0) {
      Errors::Message msg;
      msg << "Flow PK: water retention models do not cover the whole domain.";
      Exceptions::amanzi_throw(msg);  
    }
  }
}


/* ******************************************************************
* Define a fully saturated soil. 
****************************************************************** */
void RelativePermeability::SetFullySaturated()
{
  Krel_cells_->PutScalar(1.0);
  Krel_faces_->PutScalar(1.0);

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Krel_amanzi_.clear();
  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    std::vector<double> krel(nfaces, 1.0);
    Krel_amanzi_.push_back(krel);
  }
}

}  // namespace AmanziFlow
}  // namespace Amanzi

