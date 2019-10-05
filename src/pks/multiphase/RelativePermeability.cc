/*
This is the multiphase component of the Amanzi code. 

Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
         Quan Bui (mquanbui@math.umd.edu)
This file implements the main methods of the class for RelativePermeability. 
Unlike the one in Flow, this is specialized for multiphase. The relative 
permeability is computed from water saturation, instead of pressure (as in Richards).
*/

#include "errors.hh"

#include "MultiphaseDefs.hh"
#include "RelativePermeability.hh"
#include "WaterRetentionModel.hh"
#include "WRM_Simple.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Initialize internal data.
****************************************************************** */
void RelativePermeability::Init(std::string phase_name, Teuchos::ParameterList& plist)
{
  //atm_pressure = p0;
  phase_ = phase_name;
  ProcessParameterList_(plist);

  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);
  cvs.AddComponent("face", AmanziMesh::FACE, 1);

  Krel_ = Teuchos::rcp(new CompositeVector(cvs));
  dKdS_ = Teuchos::rcp(new CompositeVector(cvs));

  //method_ = FLOW_RELATIVE_PERM_NONE;
  Krel_->PutScalarMasterAndGhosted(1.0);
  dKdS_->PutScalarMasterAndGhosted(0.0);

  pp_ = Teuchos::rcp(new ParallelCommunication(mesh_));
  map_c2mb_ = Teuchos::rcp(new Epetra_IntVector(mesh_->cell_map(true)));
  PopulateMapC2MB_();
}


/* ******************************************************************
* Defines relative permeability ONLY for cells.
****************************************************************** */
void RelativePermeability::Compute(const CompositeVector& saturation_water)
{
  const Epetra_MultiVector& Sw_cell = *saturation_water.ViewComponent("cell");
  Epetra_MultiVector& Krel_cell = *Krel_->ViewComponent("cell");
  Epetra_MultiVector& dKdS_cell = *dKdS_->ViewComponent("cell");

  for (int mb = 0; mb < WRM_.size(); mb++) {
    std::string region = WRM_[mb]->region();

    std::vector<AmanziMesh::Entity_ID> block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double Sw = Sw_cell[0][*i];
      Krel_cell[0][*i] = WRM_[mb]->k_relative(Sw, phase_);
      dKdS_cell[0][*i] = WRM_[mb]->dKdS(Sw, phase_);
      //std::cout << "cell: " << *i << "; Sw: " << Sw << "; krel: " << Krel_cell[0][*i] << "\n";
      //dKdP_cell[0][*i] = -WRM_[mb]->dKdPc(pc);  // Negative sign indicates that dKdP = -dKdPc.
    }
  }
}


/* ******************************************************************
* Use analytical formula for derivative dS/dP.                                               
****************************************************************** */
void RelativePermeability::DerivedSdP(const Epetra_MultiVector& p, Epetra_MultiVector& ds)
{
  Errors::Message msg;
  msg << "Multiphase PK: RelativePermeability::DerivedSdP is not implemented.\n";
  Exceptions::amanzi_throw(msg);
  /*
  for (int mb = 0; mb < WRM_.size(); mb++) {
    std::string region = WRM_[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double pc = atm_pressure - p[0][*i];
      ds[0][*i] = -WRM_[mb]->dSdPc(pc);  // Negative sign indicates that dSdP = -dSdPc.
    }
  } */
}


/* ******************************************************************
* Use analytical formula for derivative dK/dP.                                               
****************************************************************** */
void RelativePermeability::DerivedKdP(const Epetra_MultiVector& p, Epetra_MultiVector& dk)
{
  Errors::Message msg;
  msg << "Multiphase PK: RelativePermeability::DerivedSdP is not implemented.\n";
  Exceptions::amanzi_throw(msg);
  /*
  for (int mb = 0; mb < WRM_.size(); mb++) {
    std::string region = WRM_[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double pc = atm_pressure - p[0][*i];
      dk[0][*i] = -WRM_[mb]->dKdPc(pc);  // Negative sign indicates that dKdP = -dKdPc.
    }
  } */
}


/* ******************************************************************
* Calculates gravity flux (K * g) * n and normalizes the result.
****************************************************************** */
void RelativePermeability::ComputeGravityFlux(
    const std::vector<WhetStone::Tensor>& K, const AmanziGeometry::Point& g,
    Teuchos::RCP<CompositeVector> flux)
{
  Epetra_MultiVector& flx = *flux->ViewComponent("face", true);
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells_owned; c++) {
    AmanziGeometry::Point Kg = K[c] * g;
    Kg /= norm(Kg);
 
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    for (int n = 0; n < nfaces; n++) {
      int f = faces[n];
      const AmanziGeometry::Point& normal = mesh_->face_normal(f);
      flx[0][f] = (Kg * normal) / mesh_->face_area(f);
    }
  }

  flux->ScatterMasterToGhosted("face");
} 


/* ******************************************************************
* Auxiliary map from cells to WRM models.                                               
****************************************************************** */
void RelativePermeability::PopulateMapC2MB_()
{
  Epetra_IntVector& map = *map_c2mb_;
  map.PutValue(-1);

  for (int mb = 0; mb < WRM_.size(); mb++) {
    std::string region = WRM_[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      map[*i] = mb;
    }
  }
  
  pp_->CopyMasterCell2GhostCell(map);

  // internal check
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int c = 0; c < ncells_owned; c++) {
    if (map[c] < 0) {
      Errors::Message msg;
      msg << "Multiphase PK: water retention models do not cover the whole domain.\n";
      Exceptions::amanzi_throw(msg);  
    }
  }
}

}  // namespace Flow
}  // namespace Amanzi

