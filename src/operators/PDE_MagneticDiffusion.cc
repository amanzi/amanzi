/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The implementation assumes that underlying discretization is 
  of serendipity type, has dofs on primary entities.
*/

#include <vector>

// TPLs
#include "Epetra_Vector.h"

// Amanzi
#include "errors.hh"
#include "MFD3D_Electromagnetics.hh"
#include "VEM_NedelecSerendipityType2.hh"

// Amanzi::Operators
#include "PDE_MagneticDiffusion.hh"
#include "Op.hh"

namespace Amanzi {
namespace Operators {

/* ******************************************************************
* Populate contains of elemental matrices.
* NOTE: The input parameters are not yet used.
****************************************************************** */
void PDE_MagneticDiffusion::UpdateMatrices(
    const Teuchos::Ptr<const CompositeVector>& u,
    const Teuchos::Ptr<const CompositeVector>& p)
{
  WhetStone::DenseMatrix Mcell, Ccell, Acell;

  WhetStone::Tensor Kc(mesh_->space_dimension(), 1);
  Kc(0, 0) = 1.0;
  
  for (int c = 0; c < ncells_owned; c++) {
    auto tmp = Teuchos::rcp_dynamic_cast<WhetStone::MFD3D_Electromagnetics>(mfd_);
    if (tmp != Teuchos::null) {
      tmp->StiffnessMatrix(c, Kc, Acell, Mcell, Ccell);
    } else {
      auto tmp2 = Teuchos::rcp_dynamic_cast<WhetStone::VEM_NedelecSerendipityType2>(mfd_);
      tmp2->StiffnessMatrix(c, Kc, Acell, Mcell, Ccell);
    }

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
  B.ScatterMasterToGhosted("face");
  global_op_->rhs()->PutScalarGhosted(0.0);

  const Epetra_MultiVector& Bf = *B.ViewComponent("face", true);
  Epetra_MultiVector& rhs_e = *global_op_->rhs()->ViewComponent("edge", true);

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
    mfd_->MassMatrix(c, (*K_)[c], Tcell);
    Acell += Tcell;

    int ndofs_e = Tcell.NumRows();
    int ndofs_f = Mcell.NumRows();
    WhetStone::DenseVector v1(ndofs_f), v2(ndofs_f), v3(ndofs_e);

    int nde = ndofs_e / nedges;  // assume serendipity method
    int ndf = ndofs_f / nfaces;

    int m(0);
    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      for (int k = 0; k < ndf; ++k) {
        v1(m++) = Bf[k][f] * dirs[n] * mesh_->face_area(f);
      }
    }

    Mcell.Multiply(v1, v2, false);
    Ccell.Multiply(v2, v3, true);

    m = 0;
    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      for (int k = 0; k < nde; ++k) {
        rhs_e[k][e] += v3(m++);
      }
    }
  }

  global_op_->rhs()->GatherGhostedToMaster("edge", Add);
}


/* ******************************************************************
* Solution postprocessing
* **************************************************************** */
void PDE_MagneticDiffusion::ModifyFields(
   CompositeVector& E, CompositeVector& B, double dt)
{
  E.ScatterMasterToGhosted("edge");

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

    int ndofs_e = Ccell.NumCols();
    int ndofs_f = Ccell.NumRows();
    int nde = ndofs_e / nedges;  // assumes serendipity method
    int ndf = ndofs_f / nfaces;

    WhetStone::DenseVector v1(ndofs_e), v2(ndofs_f);

    int m(0);
    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      for (int k = 0; k < nde; ++k) {
        v1(m++) = Ee[k][e];
      }
    }

    Ccell.Multiply(v1, v2, false);

    m = 0;
    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      if (!fflag[f]) {
        for (int k = 0; k < ndf; ++k) {
          Bf[k][f] -= dt * v2(m++) * dirs[n] / mesh_->face_area(f);
        }
        fflag[f] = true;
      } else {
        m += ndf;
      }
    }
  }
}


/* ******************************************************************
* Calculates Ohmic heating
****************************************************************** */
double PDE_MagneticDiffusion::CalculateOhmicHeating(const CompositeVector& E)
{
  const Epetra_MultiVector& Ee = *E.ViewComponent("edge", true);
  E.ScatterMasterToGhosted("edge");

  double energy(0.0);
  AmanziMesh::Entity_ID_List edges;

  for (int c = 0; c < ncells_owned; ++c) {
    mesh_->cell_get_edges(c, &edges);
    int nedges = edges.size();

    WhetStone::DenseMatrix Mcell;
    mfd_->MassMatrix(c, (*K_)[c], Mcell);

    int ndofs_e = Mcell.NumRows();
    int nde = ndofs_e / nedges;
    WhetStone::DenseVector v1(ndofs_e), v2(ndofs_e);

    int m(0);
    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];
      for (int k = 0; k < nde; ++k) {
        v1(m++) = Ee[k][e];
      }
    }

    Mcell.Multiply(v1, v2, false);
    energy += v1 * v2;
  }

  // parallel collective operation
  double tmp(energy);
  mesh_->get_comm()->SumAll(&tmp, &energy, 1);

  return energy;
}


/* ******************************************************************
* Calculates integral of 1/2 |B|^2
****************************************************************** */
double PDE_MagneticDiffusion::CalculateMagneticEnergy(const CompositeVector& B)
{
  const Epetra_MultiVector& Bf = *B.ViewComponent("face", true);
  B.ScatterMasterToGhosted("face");

  std::vector<int> dirs;
  AmanziMesh::Entity_ID_List faces;

  double energy(0.0);
  for (int c = 0; c < ncells_owned; ++c) {
    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    int nfaces = faces.size();

    const WhetStone::DenseMatrix& Mcell = mass_op_[c];
    int ndofs_f = Mcell.NumRows();
    int ndf = ndofs_f / nfaces;  // serendipity assumption

    WhetStone::DenseVector v1(ndofs_f), v2(ndofs_f);

    int m(0);
    for (int n = 0; n < nfaces; ++n) {
      int f = faces[n];
      for (int k = 0; k < ndf; ++k) {
        v1(m++) = Bf[k][f] * dirs[n] * mesh_->face_area(f);
      }
    }

    Mcell.Multiply(v1, v2, false);
    energy += v1 * v2;
  }

  // parallel collective operation
  double tmp(energy);
  mesh_->get_comm()->SumAll(&tmp, &energy, 1);

  return energy / 2;
}


/* ******************************************************************
* Useful tools
* **************************************************************** */
double PDE_MagneticDiffusion::CalculateDivergence(
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
void PDE_MagneticDiffusion::InitMagneticDiffusion_(Teuchos::ParameterList& plist)
{ 
  mass_op_.resize(ncells_owned);
  curl_op_.resize(ncells_owned);
}

}  // namespace Operators
}  // namespace Amanzi
