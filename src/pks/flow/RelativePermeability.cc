/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "errors.hh"

#include "FlowDefs.hh"
#include "RelativePermeability.hh"
#include "WRM.hh"
#include "WRM_vanGenuchten.hh"
#include "WRM_BrooksCorey.hh"
#include "WRM_fake.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Initialize internal data.
****************************************************************** */
void RelativePermeability::Init(double p0, Teuchos::ParameterList& plist)
{
  atm_pressure = p0;
  ProcessParameterList_(plist);

  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);
  cvs.AddComponent("face", AmanziMesh::FACE, 1);

  Krel_ = Teuchos::rcp(new CompositeVector(cvs));
  dKdP_ = Teuchos::rcp(new CompositeVector(cvs));

  method_ = FLOW_RELATIVE_PERM_NONE;
  Krel_->PutScalarMasterAndGhosted(1.0);
  dKdP_->PutScalarMasterAndGhosted(0.0);

  pp_ = Teuchos::rcp(new ParallelCommunication(mesh_));
  map_c2mb_ = Teuchos::rcp(new Epetra_IntVector(mesh_->cell_map(true)));
  PopulateMapC2MB_();
}


/* ******************************************************************
* Defines relative permeability ONLY for cells.
****************************************************************** */
void RelativePermeability::Compute(const CompositeVector& pressure)
{
  const Epetra_MultiVector& p = *pressure.ViewComponent("cell");
  Epetra_MultiVector& Krel_cell = *Krel_->ViewComponent("cell");
  Epetra_MultiVector& dKdP_cell = *dKdP_->ViewComponent("cell");

  for (int mb = 0; mb < wrm_.size(); mb++) {
    std::string region = wrm_[mb]->region();

    std::vector<AmanziMesh::Entity_ID> block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double pc = atm_pressure - p[0][*i];
      Krel_cell[0][*i] = wrm_[mb]->k_relative(pc);
      dKdP_cell[0][*i] = -wrm_[mb]->dKdPc(pc);  // Negative sign indicates that dKdP = -dKdPc.
    }
  }
}


/* ******************************************************************
* Use analytical formula for derivative dS/dP.                                               
****************************************************************** */
void RelativePermeability::DerivedSdP(const Epetra_MultiVector& p, Epetra_MultiVector& ds)
{
  for (int mb = 0; mb < wrm_.size(); mb++) {
    std::string region = wrm_[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double pc = atm_pressure - p[0][*i];
      ds[0][*i] = -wrm_[mb]->dSdPc(pc);  // Negative sign indicates that dSdP = -dSdPc.
    }
  }
}


/* ******************************************************************
* Use analytical formula for derivative dK/dP.                                               
****************************************************************** */
void RelativePermeability::DerivedKdP(const Epetra_MultiVector& p, Epetra_MultiVector& dk)
{
  for (int mb = 0; mb < wrm_.size(); mb++) {
    std::string region = wrm_[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) {
      double pc = atm_pressure - p[0][*i];
      dk[0][*i] = -wrm_[mb]->dKdPc(pc);  // Negative sign indicates that dKdP = -dKdPc.
    }
  }
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

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

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

  for (int mb = 0; mb < wrm_.size(); mb++) {
    std::string region = wrm_[mb]->region();
    AmanziMesh::Entity_ID_List block;
    mesh_->get_set_entities(region, AmanziMesh::CELL, AmanziMesh::OWNED, &block);

    AmanziMesh::Entity_ID_List::iterator i;
    for (i = block.begin(); i != block.end(); i++) map[*i] = mb;
  }
  
  pp_->CopyMasterCell2GhostCell(map);

  // internal check
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells_owned; c++) {
    if (map[c] < 0) {
      Errors::Message msg;
      msg << "Flow PK: water retention models do not cover the whole domain.";
      Exceptions::amanzi_throw(msg);  
    }
  }
}

}  // namespace Flow
}  // namespace Amanzi

