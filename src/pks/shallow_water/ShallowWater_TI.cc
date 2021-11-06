/*
 Shallow water PK
 
 Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
 Amanzi is released under the three-clause BSD License.
 The terms of use and "as is" disclaimer for this license are
 provided in the top-level COPYRIGHT file.
 
 Author: Svetlana Tokareva (tokareva@lanl.gov)
 */

#include "ShallowWater_PK.hh"

namespace Amanzi {
namespace ShallowWater {

void ShallowWater_PK::FunctionalTimeDerivative(double t, const TreeVector& A,
                                               TreeVector& f)
{
  bool failed = false;
  
  const Epetra_MultiVector h_temp = *A.SubVector(0)->Data()->ViewComponent("cell");
  const Epetra_MultiVector q_temp = *A.SubVector(1)->Data()->ViewComponent("cell");
  
  Epetra_MultiVector& h_new = *f.SubVector(0)->Data()->ViewComponent("cell");
  Epetra_MultiVector& q_new = *f.SubVector(1)->Data()->ViewComponent("cell");
  
  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  
  double tmp[1];
  S_->GetConstantVectorData("gravity", "state")->Norm2(tmp);
  double g = tmp[0];

  S_->GetFieldEvaluator(discharge_key_)->HasFieldChanged(S_.ptr(), passwd_);

  // distribute data to ghost cells
  S_->GetFieldData(total_depth_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(ponded_depth_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(velocity_key_)->ScatterMasterToGhosted("cell");
  S_->GetFieldData(discharge_key_)->ScatterMasterToGhosted("cell");

  // save a copy of primary and conservative fields
  Epetra_MultiVector& B_c = *S_->GetFieldData(bathymetry_key_, passwd_)->ViewComponent("cell", true);
  Epetra_MultiVector& B_n = *S_->GetFieldData(bathymetry_key_, passwd_)->ViewComponent("node", true);
  Epetra_MultiVector& h_c = *S_->GetFieldData(ponded_depth_key_, passwd_)->ViewComponent("cell", true);
  Epetra_MultiVector& ht_c = *S_->GetFieldData(total_depth_key_, passwd_)->ViewComponent("cell", true);
  Epetra_MultiVector& vel_c = *S_->GetFieldData(velocity_key_, passwd_)->ViewComponent("cell", true);
  Epetra_MultiVector& riemann_f = *S_->GetFieldData(riemann_flux_key_, passwd_)->ViewComponent("face", true);

  S_->GetFieldEvaluator(discharge_key_)->HasFieldChanged(S_.ptr(), passwd_);
  Epetra_MultiVector& q_c = *S_->GetFieldData(discharge_key_, discharge_key_)->ViewComponent("cell", true);
  
  for (int c = 0; c < ncells_owned; ++c) {
    double factor = inverse_with_tolerance(h_temp[0][c]);
    h_c[0][c] = h_temp[0][c];
    q_c[0][c] = q_temp[0][c];
    q_c[1][c] = q_temp[1][c];
    vel_c[0][c] = factor * q_temp[0][c];
    vel_c[1][c] = factor * q_temp[1][c];
    ht_c[0][c] = h_c[0][c] + B_c[0][c];
  }

  // create copies of primary fields
  Epetra_MultiVector h_c_tmp(h_c);
  Epetra_MultiVector q_c_tmp(q_c);

  h_c_tmp.PutScalar(0.0);
  q_c_tmp.PutScalar(0.0);

  // limited reconstructions
  Teuchos::ParameterList plist = sw_list_->sublist("reconstruction");
  bool use_limter = true;
  if (plist.isParameter("use limiter")) {
      use_limter = plist.get<bool>("use limiter");
  }
  if (use_limter == false) {
    auto tmp1 = S_->GetFieldData(total_depth_key_, passwd_)->ViewComponent("cell", true);
    total_depth_grad_->ComputeGradient(tmp1);
    
    auto tmp5 = S_->GetFieldData(discharge_key_, discharge_key_)->ViewComponent("cell", true);
    discharge_x_grad_->ComputeGradient(tmp5, 0);
    
    auto tmp6 = S_->GetFieldData(discharge_key_, discharge_key_)->ViewComponent("cell", true);
    discharge_y_grad_->ComputeGradient(tmp6, 1);
  }
  else {
    auto tmp1 = S_->GetFieldData(total_depth_key_, passwd_)->ViewComponent("cell", true);
    total_depth_grad_->ComputeGradient(tmp1);
    limiter_->ApplyLimiter(tmp1, 0, total_depth_grad_->gradient());
    limiter_->gradient()->ScatterMasterToGhosted("cell");
    
    auto tmp5 = S_->GetFieldData(discharge_key_, discharge_key_)->ViewComponent("cell", true);
    discharge_x_grad_->ComputeGradient(tmp5, 0);
    limiter_->ApplyLimiter(tmp5, 0, discharge_x_grad_->gradient());
    limiter_->gradient()->ScatterMasterToGhosted("cell");
    
    auto tmp6 = S_->GetFieldData(discharge_key_, discharge_key_)->ViewComponent("cell", true);
    discharge_y_grad_->ComputeGradient(tmp6, 1);
    limiter_->ApplyLimiter(tmp6, 1, discharge_y_grad_->gradient());
    limiter_->gradient()->ScatterMasterToGhosted("cell");
  }

  // update boundary conditions
  if (bcs_.size() > 0)
      bcs_[0]->Compute(t, t);

  // update source (external) terms
  for (int i = 0; i < srcs_.size(); ++i) {
    srcs_[i]->Compute(t, t);
  }

  // compute source (external) values
  // coupling submodel="rate" returns volumetric flux [m^3/s] integrated over
  // the time step in the last (the second) component of local data vector
  std::vector<double> ext_S_cell(ncells_owned, 0.0);
  for (int  i = 0; i < srcs_.size(); ++i) {
    for (auto it = srcs_[i]->begin(); it != srcs_[i]->end(); ++it) {
      int c = it->first;
      ext_S_cell[c] = it->second[0];  // data unit is [m]
    }
  }
  
  // Shallow water equations have the form
  // U_t + F_x(U) + G_y(U) = S(U)

  int dir, c1, c2;
  double h, u, v, qx, qy, factor;
  AmanziMesh::Entity_ID_List cells;

  std::vector<double> FNum_rot;  // fluxes
  std::vector<double> S;         // source term
  std::vector<double> UL(3), UR(3), U;  // local state vectors

  // Simplest flux form
  // U_i^{n+1} = U_i^n - dt/vol * (F_{i+1/2}^n - F_{i-1/2}^n) + dt * S_i

  for (int f = 0; f < nfaces_wghost; ++f) {
    double farea = mesh_->face_area(f);
    const AmanziGeometry::Point& xf = mesh_->face_centroid(f);

    mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    c1 = cells[0];
    c2 = (cells.size() == 2) ? cells[1] : -1;
    if (c1 > ncells_owned && c2 == -1) continue;
    if (c2 > ncells_owned) std::swap(c1, c2);

    AmanziGeometry::Point normal = mesh_->face_normal(f, false, c1, &dir);
    normal /= farea;

    double ht_rec = total_depth_grad_->getValue(c1, xf);
    double B_rec = BathymetryEdgeValue(f, B_n);

    if (ht_rec < B_rec) {
      ht_rec = ht_c[0][c1];
      B_rec = B_c[0][c1];
    }
    double h_rec = ht_rec - B_rec;
    failed = ErrorDiagnostics_(c1, h_rec, B_rec, ht_rec);

    double qx_rec = discharge_x_grad_->getValue(c1, xf);
    double qy_rec = discharge_y_grad_->getValue(c1, xf);

    factor = inverse_with_tolerance(h_rec);
    double vx_rec = factor * qx_rec;
    double vy_rec = factor * qy_rec;

    // rotating velocity to the face-based coordinate system
    double vn, vt;
    vn =  vx_rec * normal[0] + vy_rec * normal[1];
    vt = -vx_rec * normal[1] + vy_rec * normal[0];

    UL[0] = h_rec;
    UL[1] = h_rec * vn;
    UL[2] = h_rec * vt;

    if (c2 == -1) {
      if (bcs_.size() > 0 && bcs_[0]->bc_find(f)) {
        for (int i = 0; i < 3; ++i) {
          UR[i] = bcs_[0]->bc_value(f)[i];
        }
        vn =  UR[1] * normal[0] + UR[2] * normal[1];
        vt = -UR[1] * normal[1] + UR[2] * normal[0];
        UR[1] = UR[0] * vn;
        UR[2] = UR[0] * vt;
      }
      else {
        UR = UL;
      }

    } else {
      ht_rec = total_depth_grad_->getValue(c2, xf);

      if (ht_rec < B_rec) {
        ht_rec = ht_c[0][c2];
        B_rec = B_c[0][c2];
      }
      h_rec = ht_rec - B_rec;
      failed = ErrorDiagnostics_(c2, h_rec, B_rec, ht_rec);

      qx_rec = discharge_x_grad_->getValue(c2, xf);
      qy_rec = discharge_y_grad_->getValue(c2, xf);

      factor = inverse_with_tolerance(h_rec);
      vx_rec = factor * qx_rec;
      vy_rec = factor * qy_rec;

      vn =  vx_rec * normal[0] + vy_rec * normal[1];
      vt = -vx_rec * normal[1] + vy_rec * normal[0];

      UR[0] = h_rec;
      UR[1] = h_rec * vn;
      UR[2] = h_rec * vt;
    }

    FNum_rot = numerical_flux_->Compute(UL, UR);

    h  = FNum_rot[0];
    qx = FNum_rot[1] * normal[0] - FNum_rot[2] * normal[1];
    qy = FNum_rot[1] * normal[1] + FNum_rot[2] * normal[0];

    // save riemann mass flux
    riemann_f[0][f] = FNum_rot[0] * farea * dir;

    // add fluxes to temporary fields
    double vol = mesh_->cell_volume(c1);
    factor = farea / vol;
    h_c_tmp[0][c1] -= h * factor;
    q_c_tmp[0][c1] -= qx * factor;
    q_c_tmp[1][c1] -= qy * factor;

    if (c2 != -1) {
      vol = mesh_->cell_volume(c2);
      factor = farea / vol;
      h_c_tmp[0][c2] += h * factor;
      q_c_tmp[0][c2] += qx * factor;
      q_c_tmp[1][c2] += qy * factor;
    }
  }

  // sources (bathymetry, flux exchange, etc)
  // the code should not fail after that
  U.resize(3);

  for (int c = 0; c < ncells_owned; ++c) {
    U[0] = h_temp[0][c];
    U[1] = q_temp[0][c];
    U[2] = q_temp[1][c];
    
    S = NumericalSource(U, c);

    h  = h_c_tmp[0][c] + (S[0] + ext_S_cell[c]);
    qx = q_c_tmp[0][c] + S[1];
    qy = q_c_tmp[1][c] + S[2];

    h_new[0][c] = h;
    q_new[0][c] = qx;
    q_new[1][c] = qy;
  }
}

} // namespace ShallowWater
} // namespace Amanzi
