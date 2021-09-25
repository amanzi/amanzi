/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Relative permeability is a function of saturation: f(pc(sat)).
  Format: saturation, relative premeability, capillary pressure.
*/

#include "Mesh_Algorithms.hh"
#include "OperatorDefs.hh"

#include "FlowDefs.hh"
#include "RelPerm.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Constructor.
****************************************************************** */
RelPerm::RelPerm(Teuchos::ParameterList& plist,
                 Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                 double patm,
                 const Teuchos::RCP<WRMPartition>& wrm) :
    mesh_(mesh),
    patm_(patm),
    wrm_(wrm) {};


/* ******************************************************************
* Compute rel perm, k=k(pc).
****************************************************************** */
void RelPerm::Compute(Teuchos::RCP<const CompositeVector> p,
                      const std::vector<int>& bc_model,
                      const std::vector<double>& bc_value,
                      const Teuchos::RCP<CompositeVector>& krel)
{
  const Epetra_MultiVector& p_c = *p->ViewComponent("cell");
  Epetra_MultiVector& krel_c = *krel->ViewComponent("cell");

  int ncells = krel_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    krel_c[0][c] = wrm_->second[(*wrm_->first)[c]]->k_relative(patm_ - p_c[0][c]);
  }

  // update Dirichlet boundary faces 
  Epetra_MultiVector& krel_f = *krel->ViewComponent("face", true);

  for (int f = 0; f != bc_model.size(); ++f) {
    if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
      int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh_, f);
      krel_f[0][f] = wrm_->second[(*wrm_->first)[c]]->k_relative(patm_ - bc_value[f]);
    }
  }
}


/* ******************************************************************
* Compute derivative of rel perm.
****************************************************************** */
void RelPerm::ComputeDerivative(Teuchos::RCP<const CompositeVector> p,
                                const std::vector<int>& bc_model,
                                const std::vector<double>& bc_value,
                                const Teuchos::RCP<CompositeVector>& dKdP)
{
  const Epetra_MultiVector& pres_c = *p->ViewComponent("cell");
  Epetra_MultiVector& derv_c = *dKdP->ViewComponent("cell");

  int ncells = derv_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    // Negative sign indicates that dKdP = -dKdPc.
    derv_c[0][c] = -wrm_->second[(*wrm_->first)[c]]->dKdPc(patm_ - pres_c[0][c]);
  }

  // add boundary face component
  Epetra_MultiVector& derv_f = *dKdP->ViewComponent("face", true);

  for (int f = 0; f != bc_model.size(); ++f) {
    if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
      int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh_, f);
      derv_f[0][f] = -wrm_->second[(*wrm_->first)[c]]->dKdPc(patm_ - bc_value[f]);
    }
  }

}


/* ******************************************************************
* Single pressure functions.
****************************************************************** */
double RelPerm::Compute(int c, double p) const {
  return wrm_->second[(*wrm_->first)[c]]->k_relative(patm_ - p);
}


double RelPerm::ComputeDerivative(int c, double p) const {
  return -wrm_->second[(*wrm_->first)[c]]->dKdPc(patm_ - p);
}


/* ******************************************************************
* Use analytical formula for derivative dS/dP.                                               
****************************************************************** */
void RelPerm::Compute_dSdP(const Epetra_MultiVector& p, Epetra_MultiVector& ds)
{
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int c = 0; c < ncells_owned; ++c) {
    // Negative sign indicates that dSdP = -dSdPc.
    double pc = patm_ - p[0][c];
    ds[0][c] = -wrm_->second[(*wrm_->first)[c]]->dSdPc(pc);
  }
}


/* ****************************************************************
* Plot water retention curves.
**************************************************************** */
void RelPerm::PlotWRMcurves()
{
  int MyPID = mesh_->cell_map(false).Comm().MyPID();
  if (MyPID != 0) return;

  int ndata(1000);
  for (int n = 0; n < wrm_->second.size(); ++n) {
    std::ofstream ofile;
    std::string filename("wrm_" + std::to_string(n) + ".txt");
    ofile.open(filename.c_str());

    double sr = wrm_->second[n]->residualSaturation();
    double ds = (1.0 - sr) / ndata;

    for (int i = 0; i < ndata; i++) {
      double sat = sr + ds * (i + 0.5);
      double pc = wrm_->second[n]->capillaryPressure(sat);
      double krel = wrm_->second[n]->k_relative(pc);
      ofile << sat << " " << krel << " " << pc << std::endl;
    }
    ofile << std::endl;
    ofile.close();
  }
}

}  // namespace Flow
}  // namespace Amanzi
