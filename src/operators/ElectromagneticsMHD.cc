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
#include "Epetra_Vector.h"

// Amanzi
#include "errors.hh"
#include "mfd3d_electromagnetics.hh"

// Amanzi::Operators
#include "ElectromagneticsMHD.hh"
#include "Op.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Calculate elemental matrices.
****************************************************************** */
void ElectromagneticsMHD::UpdateMatrices()
{
  std::vector<int> dirs;
  AmanziMesh::Entity_ID_List edges, faces;
  WhetStone::MFD3D_Electromagnetics mfd(mesh_);

  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = 1.0;
  
  for (int c = 0; c < ncells_owned; c++) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    mesh_->cell_get_edges(c, &edges);
    int nedges = edges.size();

    WhetStone::DenseMatrix Mcell(nfaces, nfaces);
    WhetStone::DenseMatrix Ccell(nfaces, nedges);
    WhetStone::DenseMatrix Acell(nedges, nedges);

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
void ElectromagneticsMHD::ModifyMatrices(
   CompositeVector& E, CompositeVector& B, double dt)
{
  const Epetra_MultiVector& Bf = *B.ViewComponent("face", true);
  Epetra_MultiVector& rhs_e = *global_op_->rhs()->ViewComponent("edge", true);

  WhetStone::MFD3D_Electromagnetics mfd(mesh_);

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
      v1(n) = Bf[0][f] * dirs[n] * mesh_->face_area(f);
    }

    Mcell.Multiply(v1, v2, false);
    Ccell.Multiply(v2, v3, true);

    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      rhs_e[0][e] += v3(n);
    }
  }
}


/* ******************************************************************
* Solution postprocessing
* **************************************************************** */
void ElectromagneticsMHD::ModifyFields(
   CompositeVector& E, CompositeVector& B, double dt)
{
  B.ScatterMasterToGhosted("face");

  Epetra_MultiVector& Ee = *E.ViewComponent("edge", true);
  Epetra_MultiVector& Bf = *B.ViewComponent("face", false);
  
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
        Bf[0][f] -= dt * v2(n) * dirs[n] / mesh_->face_area(f);
        fflag[f] = true;
      }
    }
  }
}


/* ******************************************************************
* Put here stuff that has to be done in constructor.
****************************************************************** */
double ElectromagneticsMHD::CalculateMagneticEnergy(const CompositeVector& B)
{
  const Epetra_MultiVector& Bf = *B.ViewComponent("face", true);

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
      v1(n) = Bf[0][f] * dirs[n] * mesh_->face_area(f);
    }

    Mcell.Multiply(v1, v2, false);
    energy += v1 * v2;
  }

  return energy / 2;
}


/* ******************************************************************
* Useful tools
* **************************************************************** */
double ElectromagneticsMHD::CalculateDivergence(
    int c, const CompositeVector& B)
{
  const Epetra_MultiVector& Bf = *B.ViewComponent("face", false);
  
  std::vector<int> dirs;
  AmanziMesh::Entity_ID_List faces;

  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
  int nfaces = faces.size();

  double div(0.0);
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];
    div += Bf[0][f] * dirs[n] * mesh_->face_area(f);
  }
  div /= mesh_->cell_volume(c);

  return div;
}


/* ******************************************************************
* Put here stuff that has to be done in constructor.
****************************************************************** */
void ElectromagneticsMHD::InitElectromagneticsMHD_()
{
  mass_op_.resize(ncells_owned);
  curl_op_.resize(ncells_owned);
}

}  // namespace Operators
}  // namespace Amanzi
