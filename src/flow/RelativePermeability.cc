/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "errors.hh"

#include "WaterRetentionModel.hh"
#include "WRM_vanGenuchten.hh"
#include "WRM_BrooksCorey.hh"
#include "WRM_fake.hh"

#include "RelativePermeability.hh"
#include "FlowDefs.hh"

namespace Amanzi {
namespace AmanziFlow {

/* ******************************************************************
* Initialize internal data.
****************************************************************** */
void RelativePermeability::Init(double p0, Teuchos::RCP<State> S)
{
  atm_pressure = p0;
  S_ = S;

  ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);

  nfaces_owned = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::USED);

  const Epetra_Map& cmap_wghost = mesh_->cell_map(true);
  const Epetra_Map& fmap_wghost = mesh_->face_map(true);

  const CompositeVector& pressure = *S_->GetFieldData("pressure");
  Krel_ = Teuchos::rcp(new CompositeVector(pressure));
  dKdP_ = Teuchos::rcp(new CompositeVector(pressure));

  upwind_cell = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));
  downwind_cell = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));
  face_flag = Teuchos::rcp(new Epetra_IntVector(fmap_wghost));

  method_ = FLOW_RELATIVE_PERM_NONE;
  SetFullySaturated();

  CompositeVectorSpace cvs;
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  map_c2mb_ = Teuchos::rcp(new CompositeVector(cvs));
}


/* ******************************************************************
* A wrapper for updating relative permeabilities.
****************************************************************** */
void RelativePermeability::Compute(const CompositeVector& pressure,
                                   const std::vector<int>& bc_model,
                                   const std::vector<bc_tuple>& bc_values)
{
  Epetra_MultiVector& Krel_cells = *Krel_->ViewComponent("cell", true);
  Epetra_MultiVector& Krel_faces = *Krel_->ViewComponent("face", true);

  if (method_ == FLOW_RELATIVE_PERM_UPWIND_GRAVITY ||
      method_ == FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX ||
      method_ == FLOW_RELATIVE_PERM_ARITHMETIC_MEAN) {
    ComputeOnFaces(pressure, bc_model, bc_values);
    Krel_cells.PutScalar(1.0);
    if (experimental_solver_ == FLOW_SOLVER_NEWTON || 
        experimental_solver_ == FLOW_SOLVER_PICARD_NEWTON) {
      ComputeDerivativeOnFaces(pressure, bc_model, bc_values);
    }
  } else if (method_ == FLOW_RELATIVE_PERM_AMANZI) {
    ComputeOnFaces(pressure, bc_model, bc_values);
  } else {
    ComputeInCells(pressure);
    Krel_faces.PutScalar(1.0);
  }
}

 
/* ******************************************************************
* Defines relative permeability ONLY for cells.                                               
****************************************************************** */
void RelativePermeability::ComputeInCells(const CompositeVector& pressure)
{
  const Epetra_MultiVector& p = *pressure.ViewComponent("cell");
  Epetra_MultiVector& Krel_cells = *Krel_->ViewComponent("cell");

  for (int mb = 0; mb < WRM_.size(); mb++) {
    std::string region = WRM_[mb]->region();

    std::vector<AmanziMesh::Entity_ID> block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double pc = atm_pressure - p[0][*i];
      Krel_cells[0][*i] = WRM_[mb]->k_relative(pc);
    }
  }
}


/* ******************************************************************
* Wrapper for various ways to define relative permeability of faces.
****************************************************************** */
void RelativePermeability::ComputeOnFaces(
    const CompositeVector& pressure,
    const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values)
{
  ComputeInCells(pressure);  // populates cell-based permeabilities
  Krel_->ScatterMasterToGhosted("cell");

  if (method_ == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) {  // Define K and Krel_faces
    FaceUpwindGravity_(pressure, bc_model, bc_values);

  } else if (method_ == FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX) {
    S_->GetFieldData("darcy_flux")->ScatterMasterToGhosted("face");
    const Epetra_MultiVector& flux = *S_->GetFieldData("darcy_flux")->ViewComponent("face", true);
    
    FaceUpwindFlux_(pressure, flux, bc_model, bc_values);

  } else if (method_ == FLOW_RELATIVE_PERM_ARITHMETIC_MEAN) {
    FaceArithmeticMean_(pressure);

  } else if (method_ == FLOW_RELATIVE_PERM_AMANZI) {
    FaceUpwindGravityInSoil_(pressure, bc_model, bc_values);
  }
}


/* ******************************************************************
* Defines upwinded relative permeabilities for faces using gravity. 
****************************************************************** */
void RelativePermeability::FaceUpwindGravity_(
    const CompositeVector& pressure,
    const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Epetra_MultiVector& Krel_cells = *Krel_->ViewComponent("cell", true);
  Epetra_MultiVector& Krel_faces = *Krel_->ViewComponent("face", true);
  Krel_faces.PutScalar(2.0);

  const Epetra_MultiVector& map_c2mb = *map_c2mb_->ViewComponent("cell", true);

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
          Krel_faces[0][f] = WRM_[map_c2mb[0][c]]->k_relative(pc);
        } else {
          Krel_faces[0][f] = Krel_cells[0][c];
        }
      } else {
        if (cos_angle > FLOW_RELATIVE_PERM_TOLERANCE) {
          Krel_faces[0][f] = Krel_cells[0][c]; // The upwind face.
        } else if (fabs(cos_angle) <= FLOW_RELATIVE_PERM_TOLERANCE) { 
	  if (Krel_faces[0][f] > 1.0) {
            Krel_faces[0][f] = Krel_cells[0][c] / 2;
          } else { 
            Krel_faces[0][f] += Krel_cells[0][c] / 2;  //Almost vertical face.
          }
        } 
      }
    }
  }
}


/* ******************************************************************
* Defines upwinded relative permeabilities for faces using gravity. 
****************************************************************** */
void RelativePermeability::FaceUpwindGravityInSoil_(
    const CompositeVector& pressure,
    const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Epetra_MultiVector& Krel_cells = *Krel_->ViewComponent("cell", true);
  Epetra_MultiVector& Krel_faces = *Krel_->ViewComponent("face", true);
  Krel_amanzi_.clear();

  const Epetra_MultiVector& map_c2mb = *map_c2mb_->ViewComponent("cell", true);

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
        krel[n] = (Krel_cells[0][c1] + Krel_cells[0][c2]) / 2;
        Krel_faces[0][f] = krel[n];
      } else if (flag == FLOW_PERMFLAG_INTERFACE) {
        krel[n] = Krel_cells[0][c];
        Krel_faces[0][f] = Krel_cells[0][c1];
      } else if (flag == FLOW_PERMFLAG_UPWIND) {
        if (bc_model[f] == FLOW_BC_FACE_PRESSURE && c1 < 0) {
          double pc = atm_pressure - bc_values[f][0];
          krel[n] = WRM_[map_c2mb[0][c]]->k_relative(pc);
        } else if (c1 >= 0) {
          krel[n] = Krel_cells[0][c1];
        } else {
          krel[n] = Krel_cells[0][c];
        }
        Krel_faces[0][f] = krel[n];
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
    const CompositeVector& pressure, const Epetra_MultiVector& flux,
    const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Epetra_MultiVector& Krel_cells = *Krel_->ViewComponent("cell", true);
  Epetra_MultiVector& Krel_faces = *Krel_->ViewComponent("face", true);
  Krel_faces.PutScalar(0.0);

  const Epetra_MultiVector& map_c2mb = *map_c2mb_->ViewComponent("cell", true);

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
        if (bc_model[f] == FLOW_BC_FACE_PRESSURE && flux[0][f] * dirs[n] < -tol) {
          double pc = atm_pressure - bc_values[f][0];
          Krel_faces[0][f] = WRM_[map_c2mb[0][c]]->k_relative(pc);
        } else {
          Krel_faces[0][f] = Krel_cells[0][c];   
        }
      } else {
        if (flux[0][f] * dirs[n] > tol) {
          Krel_faces[0][f] = Krel_cells[0][c];  // The upwind face.	  
        } else if (fabs(flux[0][f]) <= tol) { 
          Krel_faces[0][f] += Krel_cells[0][c] / 2;  // Almost vertical face.
        }
      }
    }
  }
}


/* ******************************************************************
* Defines relative permeabilities for faces via arithmetic averaging. 
****************************************************************** */
void RelativePermeability::FaceArithmeticMean_(const CompositeVector& pressure)
{
  AmanziMesh::Entity_ID_List cells;

  Epetra_MultiVector& Krel_cells = *Krel_->ViewComponent("cell", true);
  Epetra_MultiVector& Krel_faces = *Krel_->ViewComponent("face", true);
  Krel_faces.PutScalar(0.0);

  for (int f = 0; f < nfaces_owned; f++) {
    mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    for (int n = 0; n < ncells; n++) Krel_faces[0][f] += Krel_cells[0][cells[n]];
    Krel_faces[0][f] /= ncells;
  }
}


/* ******************************************************************
* Wrapper for various ways to define dKdP on faces.
****************************************************************** */
void RelativePermeability::ComputeDerivativeOnFaces(
    const CompositeVector& pressure,
    const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values)
{
  const Epetra_MultiVector& p = *pressure.ViewComponent("cell");
  Epetra_MultiVector& dKdP_cells = *dKdP_->ViewComponent("cell");
  Epetra_MultiVector& dKdP_faces = *dKdP_->ViewComponent("face", true);

  DerivedKdP(p, dKdP_cells);  // populates cell-based permeabilities
  dKdP_->ScatterMasterToGhosted("cell");

  if (method_ == FLOW_RELATIVE_PERM_UPWIND_GRAVITY) {
    DerivativeFaceUpwindGravity_(pressure, bc_model, bc_values);
  } else if (method_ == FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX) {
    const Epetra_MultiVector& flux = *S_->GetFieldData("darcy_flux")->ViewComponent("face", true);
    DerivativeFaceUpwindFlux_(pressure, flux, bc_model, bc_values);
  }
}


/* ******************************************************************
* Defines upwind value of dKdP on faces using gravity. 
****************************************************************** */
void RelativePermeability::DerivativeFaceUpwindGravity_(
   const CompositeVector& pressure,
   const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Epetra_MultiVector& dKdP_cells = *dKdP_->ViewComponent("cell", true);
  Epetra_MultiVector& dKdP_faces = *dKdP_->ViewComponent("face", true);
  dKdP_faces.PutScalar(0.0);

  const Epetra_MultiVector& map_c2mb = *map_c2mb_->ViewComponent("cell", true);

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      double cos_angle = (normal * Kgravity_unit_[c]) * dirs[n] / mesh_->face_area(f);
      
      if (bc_model[f] != FLOW_BC_FACE_NULL){
        if (bc_model[f] == FLOW_BC_FACE_PRESSURE && cos_angle < -FLOW_RELATIVE_PERM_TOLERANCE) {
          int mb = map_c2mb[0][c];
          double pc = atm_pressure - bc_values[f][0];
          dKdP_faces[0][f] = WRM_[mb]->dKdPc(pc);
        } else {
          dKdP_faces[0][f] = dKdP_cells[0][c];
        }
      } else {
        if (cos_angle > FLOW_RELATIVE_PERM_TOLERANCE) {
          dKdP_faces[0][f] = dKdP_cells[0][c];  // The upwind face.
        } else if (fabs(cos_angle) <= FLOW_RELATIVE_PERM_TOLERANCE) { 
          dKdP_faces[0][f] += dKdP_cells[0][c] / 2;  // Almost vertical face.
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
    const CompositeVector& pressure, const Epetra_MultiVector& flux,
    const std::vector<int>& bc_model, const std::vector<bc_tuple>& bc_values)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  Epetra_MultiVector& dKdP_cells = *dKdP_->ViewComponent("cell", true);
  Epetra_MultiVector& dKdP_faces = *dKdP_->ViewComponent("face", true);
  dKdP_faces.PutScalar(0.0);

  const Epetra_MultiVector& map_c2mb = *map_c2mb_->ViewComponent("cell", true);

  double max_flux;
  flux.MaxValue(&max_flux);
  double tol = FLOW_RELATIVE_PERM_TOLERANCE * max_flux;

  for (int c = 0; c < ncells_wghost; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];

      if (bc_model[f] != FLOW_BC_FACE_NULL) {  // The boundary face.
        if (bc_model[f] == FLOW_BC_FACE_PRESSURE && flux[0][f] * dirs[n] < -tol) {
          int mb = map_c2mb[0][c];
          double pc = atm_pressure - bc_values[f][0];
          dKdP_faces[0][f] = WRM_[mb]->dKdPc(pc);
        } else {
          dKdP_faces[0][f] = dKdP_cells[0][c];
        }
      } else {
        if (flux[0][f] * dirs[n] > tol) {
          dKdP_faces[0][f] = dKdP_cells[0][c];  // The upwind face.
        } else if (fabs(flux[0][f]) <= tol) { 
          dKdP_faces[0][f] += dKdP_cells[0][c] / 2; // Zero flux face.
        }
      }
    }
  }
}


/* ******************************************************************
* Use analytical formula for derivative dS/dP.                                               
****************************************************************** */
void RelativePermeability::DerivedSdP(const Epetra_MultiVector& p, Epetra_MultiVector& ds)
{
  for (int mb = 0; mb < WRM_.size(); mb++) {
    std::string region = WRM_[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double pc = atm_pressure - p[0][*i];
      ds[0][*i] = -WRM_[mb]->dSdPc(pc);  // Negative sign indicates that dSdP = -dSdPc.
    }
  }
}


/* ******************************************************************
* Use analytical formula for derivative dK/dP.                                               
****************************************************************** */
void RelativePermeability::DerivedKdP(const Epetra_MultiVector& p, Epetra_MultiVector& dk)
{
  for (int mb = 0; mb < WRM_.size(); mb++) {
    std::string region = WRM_[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double pc = atm_pressure - p[0][*i];
      dk[0][*i] = -WRM_[mb]->dKdPc(pc);  // Negative sign indicates that dKdP = -dKdPc.
    }
  }
}


/* ******************************************************************
* Defines upwinded relative permeabilities for faces using gravity.
* The unwind is restricted to each soil. 
****************************************************************** */
void RelativePermeability::FaceUpwindGravityInit_()
{
  const Epetra_MultiVector& map_c2mb = *map_c2mb_->ViewComponent("cell", true);

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
      int mb1 = map_c2mb[0][c1];
      int mb2 = map_c2mb[0][c2];
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
  int dim = g.dim();
  CompositeVectorSpace cvs;

  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, dim);

  CompositeVector Kg_cv(cvs);  // temporary vector
  Epetra_MultiVector& Kg_data = *Kg_cv.ViewComponent("cell");

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells_owned; c++) {
    AmanziGeometry::Point Kg = K[c] * g;
    double Kg_norm = norm(Kg);
 
    for (int i = 0; i < dim; i++) Kg_data[i][c] = Kg[i] / Kg_norm;
  }

  Kg_cv.ScatterMasterToGhosted("cell");

  int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
  Kgravity_unit_.clear();

  for (int c = 0; c < ncells_wghost; c++) {
    AmanziGeometry::Point Kg(dim); 
    for (int i = 0; i < dim; i++) Kg[i] = Kg_data[i][c];
    Kgravity_unit_.push_back(Kg);
  }

  FaceUpwindGravityInit_();
} 


/* ******************************************************************
* Auxiliary map from cells to WRM models.                                               
****************************************************************** */
void RelativePermeability::PopulateMapC2MB()
{
  Epetra_MultiVector& map_c2mb = *map_c2mb_->ViewComponent("cell", true);
  map_c2mb.PutScalar(-1);

  for (int mb = 0; mb < WRM_.size(); mb++) {
    std::string region = WRM_[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) map_c2mb[0][*i] = mb;
  }
  
  map_c2mb_->ScatterMasterToGhosted("cell");

  // internal check
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells_owned; c++) {
    if (map_c2mb[0][c] < 0) {
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
  Krel_->ViewComponent("cell", true)->PutScalar(1.0);
  Krel_->ViewComponent("face", true)->PutScalar(1.0);

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

