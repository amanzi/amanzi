/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <vector>

// TPLs

// Amanzi
#include "errors.hh"
#include "MFD3D_Electromagnetics.hh"

// Amanzi::Operators
#include "PDE_MagneticDiffusion.hh"
#include "Op.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Populate contains of elemental matrices.
****************************************************************** */
void PDE_MagneticDiffusion::UpdateMatrices()
{
  Teuchos::ParameterList plist;
  WhetStone::MFD3D_Electromagnetics mfd(plist, mesh_);
  WhetStone::DenseMatrix Mcell, Ccell, Acell;

  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = 1.0;
  
  for (int c = 0; c < ncells_owned; c++) {
    mfd.StiffnessMatrix(c, Kc, Acell, Mcell, Ccell);

    local_op_->matrices[c] = Acell;
    mass_op_[c] = Mcell;
    curl_op_[c] = Ccell;
  }
}


/* ******************************************************************
* System modification before solving the problem:
* A := invK I + dt/2 A   and  f += Curl M B 
* **************************************************************** */
void PDE_MagneticDiffusion::ModifyMatrices(
   CompositeVector& E, CompositeVector& B, double dt)
{
  B.scatterMasterToGhosted("face");
  global_op_->rhs()->PutScalarGhosted(0.0);

  const Epetra_MultiVector& Bf = *B.viewComponent("face", true);
  Epetra_MultiVector& rhs_e = *global_op_->rhs()->viewComponent("edge", true);

  Teuchos::ParameterList plist;
  WhetStone::MFD3D_Electromagnetics mfd(plist, mesh_);

  std::vector<int> dirs;
  AmanziMesh::Entity_ID_List faces, edges;

  for (int c = 0; c < ncells_owned; ++c) {
    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];
    Acell.Scale(dt / 2);

    const WhetStone::DenseMatrix& Mcell = mass_op_[c];
    const WhetStone::DenseMatrix& Ccell = curl_op_[c];

    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    mesh_->cell_get_edges(c, &edges);

    int nfaces = faces.size();
    int nedges = edges.size();

    WhetStone::DenseMatrix Tcell(nedges, nedges);
    mfd.MassMatrix(c, (*K_)[c], Tcell);
    Acell += Tcell;

    WhetStone::DenseVector v1(nfaces), v2(nfaces), v3(nedges);

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      v1(n) = Bf[0][f] * dirs[n] * mesh_->getFaceArea(f);
    }

    Mcell.Multiply(v1, v2, false);
    Ccell.Multiply(v2, v3, true);

    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      rhs_e[0][e] += v3(n);
    }
  }

  global_op_->rhs()->gatherGhostedToMaster("edge", Add);
}


/* ******************************************************************
* Solution postprocessing
* **************************************************************** */
void PDE_MagneticDiffusion::ModifyFields(
   CompositeVector& E, CompositeVector& B, double dt)
{
  E.scatterMasterToGhosted("edge");

  Epetra_MultiVector& Ee = *E.viewComponent("edge", true);
  Epetra_MultiVector& Bf = *B.viewComponent("face", false);
  
  std::vector<int> dirs;
  AmanziMesh::Entity_ID_List faces, edges;

  std::vector<bool> fflag(nedges_wghost, false);

  for (int c = 0; c < ncells_owned; ++c) {
    const WhetStone::DenseMatrix& Ccell = curl_op_[c];

    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    mesh_->cell_get_edges(c, &edges);

    int nfaces = faces.size();
    int nedges = edges.size();

    WhetStone::DenseVector v1(nedges), v2(nfaces);

    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      v1(n) = Ee[0][e];
    }

    Ccell.Multiply(v1, v2, false);

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      if (!fflag[f]) {
        Bf[0][f] -= dt * v2(n) * dirs[n] / mesh_->getFaceArea(f);
        fflag[f] = true;
      }
    }
  }
}


/* ******************************************************************
* Calculates Ohmic heating
****************************************************************** */
double PDE_MagneticDiffusion::CalculateOhmicHeating(const CompositeVector& E)
{
  const Epetra_MultiVector& Ee = *E.viewComponent("edge", true);
  E.scatterMasterToGhosted("edge");

  AmanziMesh::Entity_ID_List edges;
  Teuchos::ParameterList plist;
  WhetStone::MFD3D_Electromagnetics mfd(plist, mesh_);

  double energy(0.0);
  for (int c = 0; c < ncells_owned; ++c) {
    mesh_->cell_get_edges(c, &edges);
    int nedges = edges.size();

    WhetStone::DenseMatrix Mcell;
    mfd.MassMatrix(c, (*K_)[c], Mcell);

    WhetStone::DenseVector v1(nedges), v2(nedges);

    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      v1(n) = Ee[0][e];
    }

    Mcell.Multiply(v1, v2, false);
    energy += v1 * v2;
  }

  // parallel collective operation
  double tmp(energy);
  mesh_->getComm()->SumAll(&tmp, &energy, 1);

  return energy;
}


/* ******************************************************************
* Calculates integral of 1/2 |B|^2
****************************************************************** */
double PDE_MagneticDiffusion::CalculateMagneticEnergy(const CompositeVector& B)
{
  const Epetra_MultiVector& Bf = *B.viewComponent("face", true);
  B.scatterMasterToGhosted("face");

  std::vector<int> dirs;
  AmanziMesh::Entity_ID_List faces;

  double energy(0.0);
  for (int c = 0; c < ncells_owned; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    const WhetStone::DenseMatrix& Mcell = mass_op_[c];

    WhetStone::DenseVector v1(nfaces), v2(nfaces);

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      v1(n) = Bf[0][f] * dirs[n] * mesh_->getFaceArea(f);
    }

    Mcell.Multiply(v1, v2, false);
    energy += v1 * v2;
  }

  // parallel collective operation
  double tmp(energy);
  mesh_->getComm()->SumAll(&tmp, &energy, 1);

  return energy / 2;
}


/* ******************************************************************
* Useful tools
* **************************************************************** */
double PDE_MagneticDiffusion::CalculateDivergence(
    int c, const CompositeVector& B)
{
  const Epetra_MultiVector& Bf = *B.viewComponent("face", false);
  
  std::vector<int> dirs;
  AmanziMesh::Entity_ID_List faces;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  double div(0.0);
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    div += Bf[0][f] * dirs[n] * mesh_->getFaceArea(f);
  }
  div /= mesh_->getCellVolume(c);

  return div;
}


/* ******************************************************************
* Put here stuff that has to be done in constructor.
****************************************************************** */
void PDE_MagneticDiffusion::InitMagneticDiffusion_(Teuchos::ParameterList& plist)
{
  // Primary discretization methods
  std::string primary = plist.get<std::string>("discretization primary");

  if (primary == "mfd: default") {
    mfd_primary_ = WhetStone::ELECTROMAGNETICS_DEFAULT;
  } else if (primary == "mfd: generalized") {
    mfd_primary_ = WhetStone::ELECTROMAGNETICS_GENERALIZED;
  } else {
    Errors::Message msg;
    msg << "Primary discretization method \"" << primary << "\" is not supported.";
    Exceptions::amanzi_throw(msg);
  } 

  // Define stencil for the MFD diffusion method.
  mass_op_.resize(ncells_owned);
  curl_op_.resize(ncells_owned);
}

}  // namespace Operators
}  // namespace Amanzi
