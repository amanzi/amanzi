/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Relative permeability is a function of saturation: f(pc(sat)).
*/

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
    wrm_(wrm)
{
  Init_(plist);
};


/* ******************************************************************
* Compute rel perm, k=k(pc).
****************************************************************** */
void RelPerm::Compute(Teuchos::RCP<const CompositeVector> p,
                      Teuchos::RCP<CompositeVector> krel)
{
  const Epetra_MultiVector& p_c = *p->ViewComponent("cell", false);
  Epetra_MultiVector& krel_c = *krel->ViewComponent("cell", false);

  int ncells = krel_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    krel_c[0][c] = wrm_->second[(*wrm_->first)[c]]->k_relative(patm_ - p_c[0][c]);
  }
}


/* ******************************************************************
* Compute derivative of rel perm.
****************************************************************** */
void RelPerm::ComputeDerivative(Teuchos::RCP<const CompositeVector> p,
                                Teuchos::RCP<CompositeVector> dKdP)
{
  const Epetra_MultiVector& pres_c = *p->ViewComponent("cell", false);
  Epetra_MultiVector& derv_c = *dKdP->ViewComponent("cell", false);

  int ncells = derv_c.MyLength();
  for (int c = 0; c != ncells; ++c) {
    // Negative sign indicates that dKdP = -dKdPc.
    derv_c[0][c] = -wrm_->second[(*wrm_->first)[c]]->dKdPc(patm_ - pres_c[0][c]);
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
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (int c = 0; c < ncells_owned; ++c) {
    // Negative sign indicates that dSdP = -dSdPc.
    double pc = patm_ - p[0][c];
    ds[0][c] = -wrm_->second[(*wrm_->first)[c]]->dSdPc(pc);
  }
}


/* ******************************************************************
* Calculates gravity flux (K * g) * n and normalizes the result.
****************************************************************** */
void RelPerm::ComputeGravityFlux(
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


/* ****************************************************************
* Process string for the relative permeability
**************************************************************** */
void RelPerm::Init_(Teuchos::ParameterList& plist)
{
  std::string name = plist.get<std::string>("relative permeability");

  if (name == "upwind-gravity") {
    method_ = Flow::FLOW_RELATIVE_PERM_UPWIND_GRAVITY;
  } else if (name == "upwind-darcy_velocity") {
    method_ = Flow::FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX;
  } else if (name == "upwind-amanzi") {
    method_ = Flow::FLOW_RELATIVE_PERM_AMANZI_MFD;
  } else if (name == "other-arithmetic_average") {
    method_ = Flow::FLOW_RELATIVE_PERM_ARITHMETIC_AVERAGE;
  } else if (name == "other-harmonic_average") {
    method_ = Flow::FLOW_RELATIVE_PERM_HARMONIC_AVERAGE;
  } else {
    Errors::Message msg;
    msg << "Flow PK: unknown relative permeability method has been specified." ;
    Exceptions::amanzi_throw(msg);
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
  for (int i = 0; i < wrm_->second.size(); ++i) {
    std::stringstream fname;
    fname << "wrm_curves_" << i << ".txt";
    std::ofstream ofile;
    ofile.open(fname.str().c_str());

    double sr = wrm_->second[i]->residualSaturation();
    double ds = (1.0 - sr) / ndata;

    for (int i = 0; i < ndata; i++) {
      double sat = sr + ds * (i + 0.5);
      double pc = wrm_->second[i]->capillaryPressure(sat);
      double krel = wrm_->second[i]->k_relative(pc);
      ofile << sat << " " << krel << " " << pc << std::endl;
    }
    ofile << std::endl;
    ofile.close();
  }
}

}  // namespace Flow
}  // namespace Amanzi
