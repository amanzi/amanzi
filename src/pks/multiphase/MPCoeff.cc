/*
  MultiPhase

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Quan Bui (mquanbui@math.umd.edu)

  This file implements the main methods of the class for MPCoeff. 
  Unlike the one in Flow, this is specialized for multiphase. The relative 
  permeability is computed from water saturation, instead of pressure (as in Richards).
*/

#include "errors.hh"

#include "ModelMeshPartition.hh"
#include "MultiphaseDefs.hh"
#include "MPCoeff.hh"

namespace Amanzi {
namespace Multiphase {

/* ******************************************************************
* Initialize internal data.
****************************************************************** */
void MPCoeff::Init(std::string phase_name, Teuchos::ParameterList& plist)
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
  cvs.AddComponent("dirichlet_faces", AmanziMesh::BOUNDARY_FACE, 1);

  Krel_ = Teuchos::rcp(new CompositeVector(cvs));
  dKdS_ = Teuchos::rcp(new CompositeVector(cvs));
  rho_ = Teuchos::rcp(new CompositeVector(cvs));
  mpCoeff_ = Teuchos::rcp(new CompositeVector(cvs));
  rhoDerivKrel_ = Teuchos::rcp(new CompositeVector(cvs));
  dPc_ = Teuchos::rcp(new CompositeVector(cvs));

  Krel_->PutScalarMasterAndGhosted(1.0);
  dKdS_->PutScalarMasterAndGhosted(0.0);
  rho_->PutScalarMasterAndGhosted(1.0);
  mpCoeff_->PutScalarMasterAndGhosted(1.0);
  rhoDerivKrel_->PutScalarMasterAndGhosted(1.0);
  dPc_->PutScalarMasterAndGhosted(1.0);

  Cg_ = 1.0;
}


void MPCoeff::Init(std::string phase_name, Teuchos::ParameterList& plist, double Cg)
{
  Init(phase_name, plist);
  Cg_ = Cg;
}


/* ******************************************************************
* Defines relative permeability ONLY for cells.
****************************************************************** */
void MPCoeff::Compute(
    const CompositeVector& Sw,
    const std::vector<int>& bc_model, const std::vector<double>& bc_value)
{
  const Epetra_MultiVector& Sw_cell = *Sw.ViewComponent("cell");
  Epetra_MultiVector& Krel_cell = *Krel_->ViewComponent("cell");
  Epetra_MultiVector& dKdS_cell = *dKdS_->ViewComponent("cell");

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int c = 0; c != ncells; ++c) {
    int mb = (*wrm_->first)[c];
    double Sw = Sw_cell[0][c];
    Krel_cell[0][c] = (wrm_->second)[mb]->k_relative(Sw, phase_);
    dKdS_cell[0][c] = (wrm_->second)[mb]->dKdS(Sw, phase_);
  }

  // add boundary face component
  Krel_->ViewComponent("dirichlet_faces", true)->PutScalar(0.0);
  dKdS_->ViewComponent("dirichlet_faces", true)->PutScalar(0.0);
}


/* ******************************************************************
* Defines relative permeability ONLY for cells.
****************************************************************** */
void MPCoeff::Compute(
    const CompositeVector& primary_var, const CompositeVector& Sw,
    const std::vector<int>& bc_model, const std::vector<double>& bc_value)
{
  const Epetra_MultiVector& Sw_cell = *Sw.ViewComponent("cell");
  const Epetra_MultiVector& primary_var_cell = *primary_var.ViewComponent("cell");
  Epetra_MultiVector& Krel_cell = *Krel_->ViewComponent("cell");
  Epetra_MultiVector& dKdS_cell = *dKdS_->ViewComponent("cell");
  Epetra_MultiVector& rho_cell = *rho_->ViewComponent("cell");
  Epetra_MultiVector& mpCoeff_cell = *mpCoeff_->ViewComponent("cell");
  Epetra_MultiVector& rhoDerivKrel_cell = *rhoDerivKrel_->ViewComponent("cell");
  Epetra_MultiVector& dPc_dS_cell = *dPc_->ViewComponent("cell");

  int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int c = 0; c != ncells; ++c) {
    int mb = (*wrm_->first)[c];
    double Sw = Sw_cell[0][c];

    double tmpRhoCell;
    if (phase_ == "gas")
      tmpRhoCell = Cg_ * (primary_var_cell[0][c] + (wrm_->second)[mb]->capillaryPressure(Sw));
    else
      tmpRhoCell = primary_var_cell[0][c];

    double krel = (wrm_->second)[mb]->k_relative(Sw, phase_);
    double dKdS = (wrm_->second)[mb]->dKdS(Sw, phase_);
    Krel_cell[0][c] = krel;
    dKdS_cell[0][c] = dKdS;
    rho_cell[0][c] = tmpRhoCell;
    mpCoeff_cell[0][c] = tmpRhoCell * krel;
    rhoDerivKrel_cell[0][c] = tmpRhoCell * dKdS;
    dPc_dS_cell[0][c] = (wrm_->second)[mb]->dPc_dS(Sw);
  }

  // add boundary face component
  Krel_->ViewComponent("dirichlet_faces", true)->PutScalar(0.0);
  dKdS_->ViewComponent("dirichlet_faces", true)->PutScalar(0.0);
}


void MPCoeff::ProcessParameterList_(Teuchos::ParameterList& plist)
{
  Teuchos::ParameterList vlist;
  vo_ = new VerboseObject("Multiphase::Coeff", vlist); 

  auto wrm_list = Teuchos::rcp(new Teuchos::ParameterList(plist));
  wrm_ = CreateModelPartition<WRMmp>(mesh_, wrm_list, "water retention model");
}

}  // namespace Flow
}  // namespace Amanzi

