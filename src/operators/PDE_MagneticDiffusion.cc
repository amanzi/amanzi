/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

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
PDE_MagneticDiffusion::UpdateMatrices(
  const Teuchos::Ptr<const CompositeVector>& u,
  const Teuchos::Ptr<const CompositeVector>& p)
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
void
PDE_MagneticDiffusion::ModifyMatrices(CompositeVector& E, CompositeVector& B,
                                      double dt)
{
  B.ScatterMasterToGhosted("face");
  global_op_->rhs()->putScalarGhosted(0.0);

  const Epetra_MultiVector& Bf = *B.ViewComponent("face", true);
  Epetra_MultiVector& rhs_e = *global_op_->rhs()->ViewComponent("edge", true);

  Teuchos::ParameterList plist;
  WhetStone::MFD3D_Electromagnetics mfd(plist, mesh_);

  std::vector<int> dirs;
  AmanziMesh::Entity_ID_List faces, edges;

  for (int c = 0; c < ncells_owned; ++c) {
    WhetStone::DenseMatrix& Acell = local_op_->matrices[c];
    Acell.scale(dt / 2);

    const WhetStone::DenseMatrix& Mcell = mass_op_[c];
    const WhetStone::DenseMatrix& Ccell = curl_op_[c];

    mesh_->cell_get_faces_and_dirs(c, &faces, dirs);
    mesh_->cell_get_edges(c, &edges);

    int nfaces = faces.size();
    int nedges = edges.size();

    WhetStone::DenseMatrix Tcell(nedges, nedges);
    mfd.MassMatrix(c, (*K_)[c], Tcell);
    Acell += Tcell;

    WhetStone::DenseVector v1(nfaces), v2(nfaces), v3(nedges);

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      v1(n) = Bf[0][f] * dirs[n] * mesh_->face_area(f);
    }

    Mcell.elementWiseMultiply(v1, v2, false);
    Ccell.elementWiseMultiply(v2, v3, true);

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
PDE_MagneticDiffusion::ModifyFields(CompositeVector& E, CompositeVector& B,
                                    double dt)
{
  E.ScatterMasterToGhosted("edge");

  Epetra_MultiVector& Ee = *E.ViewComponent("edge", true);
  Epetra_MultiVector& Bf = *B.ViewComponent("face", false);

  std::vector<int> dirs;
  AmanziMesh::Entity_ID_List faces, edges;

  std::vector<bool> fflag(nedges_wghost, false);

  for (int c = 0; c < ncells_owned; ++c) {
    const WhetStone::DenseMatrix& Ccell = curl_op_[c];

    mesh_->cell_get_faces_and_dirs(c, &faces, dirs);
    mesh_->cell_get_edges(c, &edges);

    int nfaces = faces.size();
    int nedges = edges.size();

    WhetStone::DenseVector v1(nedges), v2(nfaces);

    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      v1(n) = Ee[0][e];
    }

    Ccell.elementWiseMultiply(v1, v2, false);

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      if (!fflag[f]) {
        Bf[0][f] -= dt * v2(n) * dirs[n] / mesh_->face_area(f);
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

    Mcell.elementWiseMultiply(v1, v2, false);
    energy += v1 * v2;
  }

  // parallel collective operation
  double tmp(energy);
  Teuchos::reduceAll(*mesh_->get_comm(), Teuchos::REDUCE_SUM, 1, &tmp, &energy);

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

  std::vector<int> dirs;
  AmanziMesh::Entity_ID_List faces;

  double energy(0.0);
  for (int c = 0; c < ncells_owned; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, dirs);
    int nfaces = faces.size();

    const WhetStone::DenseMatrix& Mcell = mass_op_[c];

    WhetStone::DenseVector v1(nfaces), v2(nfaces);

    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      v1(n) = Bf[0][f] * dirs[n] * mesh_->face_area(f);
    }

    Mcell.elementWiseMultiply(v1, v2, false);
    energy += v1 * v2;
  }

  // parallel collective operation
  double tmp(energy);
  Teuchos::reduceAll(*mesh_->get_comm(), Teuchos::REDUCE_SUM, 1, &tmp, &energy);

  return energy / 2;
}


/* ******************************************************************
 * Useful tools
 * **************************************************************** */
double
PDE_MagneticDiffusion::CalculateDivergence(int c, const CompositeVector& B)
{
  const Epetra_MultiVector& Bf = *B.ViewComponent("face", false);

  std::vector<int> dirs;
  AmanziMesh::Entity_ID_List faces;

  mesh_->cell_get_faces_and_dirs(c, &faces, dirs);
  int nfaces = faces.size();

  double div(0.0);
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    div += Bf[0][f] * dirs[n] * mesh_->face_area(f);
  }
  div /= mesh_->cell_volume(c, false);

  return div;
}


/* ******************************************************************
 * Put here stuff that has to be done in constructor.
 ****************************************************************** */
void
PDE_MagneticDiffusion::InitMagneticDiffusion_(Teuchos::ParameterList& plist)
{
  // Primary discretization methods
  std::string primary = plist.get<std::string>("discretization primary");

  if (primary == "mfd: default") {
    mfd_primary_ = WhetStone::ELECTROMAGNETICS_DEFAULT;
  } else if (primary == "mfd: generalized") {
    mfd_primary_ = WhetStone::ELECTROMAGNETICS_GENERALIZED;
  } else {
    Errors::Message msg;
    msg << "Primary discretization method \"" << primary
        << "\" is not supported.";
    Exceptions::amanzi_throw(msg);
  }

  // Define stencil for the MFD diffusion method.
  mass_op_.resize(ncells_owned);
  curl_op_.resize(ncells_owned);
}

} // namespace Operators
} // namespace Amanzi
