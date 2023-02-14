/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

*/

#include <vector>

// TPLs
#include "Epetra_Vector.h"

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
* NOTE: The input parameters are not yet used.
****************************************************************** */
void
PDE_MagneticDiffusion::UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                                      const Teuchos::Ptr<const CompositeVector>& p)
{
  Teuchos::ParameterList plist;
  WhetStone::MFD3D_Electromagnetics mfd(plist, mesh_);
  WhetStone::DenseMatrix Mcell, Ccell, Acell;

  WhetStone::Tensor Kc(mesh_->getSpaceDimension(), 1);
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
void
PDE_MagneticDiffusion::ModifyMatrices(CompositeVector& E, CompositeVector& B, double dt)
{
  B.ScatterMasterToGhosted("face");
  global_op_->rhs()->PutScalarGhosted(0.0);

  const Epetra_MultiVector& Bf = *B.ViewComponent("face", true);
  Epetra_MultiVector& rhs_e = *global_op_->rhs()->ViewComponent("edge", true);

  Teuchos::ParameterList plist;
  WhetStone::MFD3D_Electromagnetics mfd(plist, mesh_);

  for (int c = 0; c < ncells_owned; ++c) {
    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];
    Acell.Scale(dt / 2);

    const WhetStone::DenseMatrix& Mcell = mass_op_[c];
    const WhetStone::DenseMatrix& Ccell = curl_op_[c];

    const auto& [faces,dirs] = mesh_->getCellFacesAndDirections(c);
    auto edges = mesh_->getCellEdges(c);

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

  global_op_->rhs()->GatherGhostedToMaster("edge", Add);
}


/* ******************************************************************
* Solution postprocessing
* **************************************************************** */
void
PDE_MagneticDiffusion::ModifyFields(CompositeVector& E, CompositeVector& B, double dt)
{
  E.ScatterMasterToGhosted("edge");

  Epetra_MultiVector& Ee = *E.ViewComponent("edge", true);
  Epetra_MultiVector& Bf = *B.ViewComponent("face", false);

  std::vector<bool> fflag(nedges_wghost, false);

  for (int c = 0; c < ncells_owned; ++c) {
    const WhetStone::DenseMatrix& Ccell = curl_op_[c];

    const auto& [faces,dirs] = mesh_->getCellFacesAndDirections(c);
    auto edges = mesh_->getCellEdges(c);

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
double
PDE_MagneticDiffusion::CalculateOhmicHeating(const CompositeVector& E)
{
  const Epetra_MultiVector& Ee = *E.ViewComponent("edge", true);
  E.ScatterMasterToGhosted("edge");

  Teuchos::ParameterList plist;
  WhetStone::MFD3D_Electromagnetics mfd(plist, mesh_);

  double energy(0.0);
  for (int c = 0; c < ncells_owned; ++c) {
    auto edges = mesh_->getCellEdges(c);
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
double
PDE_MagneticDiffusion::CalculateMagneticEnergy(const CompositeVector& B)
{
  const Epetra_MultiVector& Bf = *B.ViewComponent("face", true);
  B.ScatterMasterToGhosted("face");

  double energy(0.0);
  for (int c = 0; c < ncells_owned; ++c) {
    const auto& [faces,dirs] = mesh_->getCellFacesAndDirections(c);
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
double
PDE_MagneticDiffusion::CalculateDivergence(int c, const CompositeVector& B)
{
  const Epetra_MultiVector& Bf = *B.ViewComponent("face", false);

  const auto& [faces,dirs] = mesh_->getCellFacesAndDirections(c);
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
void
PDE_MagneticDiffusion::InitMagneticDiffusion_(Teuchos::ParameterList& plist)
{
  mass_op_.resize(ncells_owned);
  curl_op_.resize(ncells_owned);
}

} // namespace Operators
} // namespace Amanzi
