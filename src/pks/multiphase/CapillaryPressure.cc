/*
This is the multiphase component of the Amanzi code. 

Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
         Quan Bui (mquanbui@math.umd.edu)
This file implements the main methods of the class for CapillaryPressure. 
Unlike the one in Flow, this is specialized for multiphase. The capillary pressure 
is computed from water saturation, instead of pressure (as in Richards).
*/

#include "errors.hh"

#include "ModelMeshPartition.hh"
#include "MultiphaseDefs.hh"
#include "CapillaryPressure.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Initialize internal data.
****************************************************************** */
void CapillaryPressure::Init(Teuchos::ParameterList& plist)
{
  //atm_pressure = p0;
  ProcessParameterList_(plist);

  CompositeVectorSpace cvs; 
  cvs.SetMesh(mesh_);
  cvs.SetGhosted(true);
  cvs.SetComponent("cell", AmanziMesh::CELL, 1);
  cvs.SetOwned(false);
  cvs.AddComponent("face", AmanziMesh::FACE, 1);

  Pc_ = Teuchos::rcp(new CompositeVector(cvs));
  dPc_dS_ = Teuchos::rcp(new CompositeVector(cvs));

  Pc_->PutScalarMasterAndGhosted(0.0);
  dPc_dS_->PutScalarMasterAndGhosted(0.0);
}


/* ******************************************************************
* Defines relative permeability ONLY for cells.
****************************************************************** */
void CapillaryPressure::Compute(const CompositeVector& saturation_water)
{
  const Epetra_MultiVector& Sw_cell = *saturation_water.ViewComponent("cell");
  Epetra_MultiVector& Pc_cell = *Pc_->ViewComponent("cell");
  Epetra_MultiVector& dPc_dS_cell = *dPc_dS_->ViewComponent("cell");

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int c = 0; c != ncells; ++c) {
    int mb = (*wrm_->first)[c];
    double Sw = Sw_cell[0][c];
    Pc_cell[0][c] = (wrm_->second)[mb]->capillaryPressure(Sw);
    dPc_dS_cell[0][c] = (wrm_->second)[mb]->dPc_dS(Sw);
  }
}


/* ******************************************************************
* List to WRM models has to be provided.
****************************************************************** */
void CapillaryPressure::ProcessParameterList_(Teuchos::ParameterList& plist)
{
  Teuchos::ParameterList vlist;
  vo_ = new VerboseObject("Multiphase::CapillaryPressure", vlist); 

  auto wrm_list = Teuchos::rcp(new Teuchos::ParameterList(plist));
  wrm_ = CreateModelPartition<WRMmp>(mesh_, wrm_list, "water retention model");
}

}  // namespace Flow
}  // namespace Amanzi

