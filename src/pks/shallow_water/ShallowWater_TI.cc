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

void
ShallowWater_PK::FunctionalTimeDerivative(double t, const TreeVector& A,
                                          TreeVector& fun)
{
  bool failed = false;

  const auto& h_temp = *A.SubVector(0)->Data()->ViewComponent("cell", true);
  const auto& q_temp = *A.SubVector(1)->Data()->ViewComponent("cell", true);

  auto& f_temp0 = *fun.SubVector(0)->Data()->ViewComponent("cell");
  auto& f_temp1 = *fun.SubVector(1)->Data()->ViewComponent("cell");

  int ncells_owned =
    mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost =
    mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost =
    mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  int nnodes_wghost =
    mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);

  // distribute data to ghost cells
  A.SubVector(0)->Data()->ScatterMasterToGhosted("cell");
  A.SubVector(1)->Data()->ScatterMasterToGhosted("cell");

  // save a copy of primary and conservative fields
  const auto& B_c =
    *S_->Get<CompositeVector>(bathymetry_key_).ViewComponent("cell", true);
  const auto& B_n =
    *S_->Get<CompositeVector>(bathymetry_key_).ViewComponent("node", true);
  auto& ht_c = *S_->GetW<CompositeVector>(total_depth_key_, passwd_)
                  .ViewComponent("cell", true);
  // auto& ht_n = *S_->GetW<CompositeVector>(total_depth_key_,
  // passwd_).ViewComponent("node", true);
  auto& vel_c = *S_->GetW<CompositeVector>(velocity_key_, passwd_)
                   .ViewComponent("cell", true);
  auto& riemann_f = *S_->GetW<CompositeVector>(riemann_flux_key_, passwd_)
                       .ViewComponent("face", true);

  for (int c = 0; c < ncells_wghost; ++c) {
    double factor = inverse_with_tolerance(h_temp[0][c]);
    vel_c[0][c] = factor * q_temp[0][c];
    vel_c[1][c] = factor * q_temp[1][c];
    ht_c[0][c] = h_temp[0][c] + B_c[0][c];
  }

  // create copies of primary fields
  Epetra_MultiVector h_c_tmp(h_temp);
  Epetra_MultiVector q_c_tmp(q_temp);

  h_c_tmp.PutScalar(0.0);
  q_c_tmp.PutScalar(0.0);

  // update boundary conditions given by [h u v]
  for (int i = 0; i < bcs_.size(); ++i) { bcs_[i]->Compute(t, t); }

  std::vector<int> bc_model(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value_hn(nnodes_wghost, 0.0);
  std::vector<double> bc_value_h(nfaces_wghost, 0.0);
  std::vector<double> bc_value_ht(nfaces_wghost, 0.0);
  std::vector<double> bc_value_qx(nfaces_wghost, 0.0);
  std::vector<double> bc_value_qy(nfaces_wghost, 0.0);

  // extract velocity and compute qx, qy, h BC at faces
  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->get_bc_name() == "ponded-depth") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int n = it->first;
        bc_value_hn[n] = it->second[0];
      }
    }
  }

  AmanziMesh::Entity_ID_List nodes;
  // extract ponded depth BC at nodes to ensure well-balancedness
  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->get_bc_name() == "velocity") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        mesh_->face_get_nodes(f, &nodes);

        bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value_h[f] = (bc_value_hn[nodes[0]] + bc_value_hn[nodes[1]]) / 2;
        bc_value_qx[f] = bc_value_h[f] * it->second[0];
        bc_value_qy[f] = bc_value_h[f] * it->second[1];
        bc_value_ht[f] =
          bc_value_h[f] + (B_n[0][nodes[0]] + B_n[0][nodes[1]]) / 2;
      }
    }
  }

  // limited reconstructions using boundary data
  auto tmp1 =
    S_->GetW<CompositeVector>(total_depth_key_, Tags::DEFAULT, passwd_)
      .ViewComponent("cell", true);
  total_depth_grad_->Compute(tmp1);
  if (use_limiter_)
    limiter_->ApplyLimiter(tmp1, 0, total_depth_grad_, bc_model, bc_value_ht);
  total_depth_grad_->data()->ScatterMasterToGhosted("cell");

  // Total depth recalculation for positivity
  //TotalDepthReconstruct();
  
  auto tmp5 = A.SubVector(1)->Data()->ViewComponent("cell", true);
  discharge_x_grad_->Compute(tmp5, 0);
  if (use_limiter_)
    limiter_->ApplyLimiter(tmp5, 0, discharge_x_grad_, bc_model, bc_value_qx);
  discharge_x_grad_->data()->ScatterMasterToGhosted("cell");

  auto tmp6 = A.SubVector(1)->Data()->ViewComponent("cell", true);
  discharge_y_grad_->Compute(tmp6, 1);
  if (use_limiter_)
    limiter_->ApplyLimiter(tmp6, 1, discharge_y_grad_, bc_model, bc_value_qy);
  discharge_y_grad_->data()->ScatterMasterToGhosted("cell");

  // update source (external) terms
  for (int i = 0; i < srcs_.size(); ++i) { srcs_[i]->Compute(t, t); }

  // compute source (external) values
  // coupling submodel="rate" returns volumetric flux [m^3/s] integrated over
  // the time step in the last (the second) component of local data vector
  std::vector<double> ext_S_cell(ncells_owned, 0.0);
  for (int i = 0; i < srcs_.size(); ++i) {
    for (auto it = srcs_[i]->begin(); it != srcs_[i]->end(); ++it) {
      int c = it->first;
      ext_S_cell[c] = it->second[0]; // data unit is [m]
    }
  }

  // Shallow water equations have the form
  // U_t + F_x(U) + G_y(U) = S(U)

  int dir, c1, c2;
  double h, qx, qy, factor;
  AmanziMesh::Entity_ID_List cells;

  std::vector<double> FNum_rot;        // fluxes
  std::vector<double> S;               // source term
  std::vector<double> UL(3), UR(3), U; // local state vectors

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
    //ht_rec = ht_cell_face[c1][f];
    double B_rec = BathymetryEdgeValue(f, B_n);
    
    const auto& xc = mesh_->cell_centroid(c1);

    if (ht_rec < B_rec && std::abs(ht_rec - B_rec) > 1.e-14) {
    	std::cout<<"time: t = "<<t<<", c = "<<c1<<", negative h 1 = : "<<ht_rec-B_rec<<" | ht_rec = "<<ht_rec<<" < "<<B_rec<<std::endl;
      ht_rec = ht_c[0][c1];
      //ht_rec = B_rec;
      B_rec = B_c[0][c1];
    } 
    double h_rec = ht_rec - B_rec;
    if (std::abs(h_rec) < 1.e-14) { h_rec = 0.0; }
    failed = ErrorDiagnostics_(t, c1, h_rec, B_rec, ht_rec);

    double qx_rec = discharge_x_grad_->getValue(c1, xf);
    double qy_rec = discharge_y_grad_->getValue(c1, xf);

    factor = inverse_with_tolerance(h_rec);
    double vx_rec = factor * qx_rec;
    double vy_rec = factor * qy_rec;

    // rotating velocity to the face-based coordinate system
    double vn, vt;
    vn = vx_rec * normal[0] + vy_rec * normal[1];
    vt = -vx_rec * normal[1] + vy_rec * normal[0];

    UL[0] = h_rec;
    UL[1] = h_rec * vn;
    UL[2] = h_rec * vt;

    if (c2 == -1) {
      if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
        UR[0] = bc_value_h[f];
        UR[1] = bc_value_qx[f] * normal[0] + bc_value_qy[f] * normal[1];
        UR[2] = -bc_value_qx[f] * normal[1] + bc_value_qy[f] * normal[0];
      } else {
        // default outflow BC
        UR = UL;
      }
    } else {
      ht_rec = total_depth_grad_->getValue(c2, xf);
      //ht_rec = ht_cell_face[c2][f];

      if (ht_rec < B_rec && std::abs(ht_rec - B_rec) > 1.e-14) {
      	std::cout<<"time: t = "<<t<<", c = "<<c2<<", negative h 2 = : "<<ht_rec-B_rec<<" | ht_rec = "<<ht_rec<<" < "<<B_rec<<std::endl;
        ht_rec = ht_c[0][c2];
        //ht_rec = B_rec;
        B_rec = B_c[0][c2];
      }
      h_rec = ht_rec - B_rec;
      if (std::abs(h_rec) < 1.e-14) { h_rec = 0.0; }
      failed = ErrorDiagnostics_(t, c2, h_rec, B_rec, ht_rec);

      qx_rec = discharge_x_grad_->getValue(c2, xf);
      qy_rec = discharge_y_grad_->getValue(c2, xf);

      factor = inverse_with_tolerance(h_rec);
      vx_rec = factor * qx_rec;
      vy_rec = factor * qy_rec;

      vn = vx_rec * normal[0] + vy_rec * normal[1];
      vt = -vx_rec * normal[1] + vy_rec * normal[0];

      UR[0] = h_rec;
      UR[1] = h_rec * vn;
      UR[2] = h_rec * vt;
    }

    FNum_rot = numerical_flux_->Compute(UL, UR);

    h = FNum_rot[0];
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

    h = h_c_tmp[0][c] + (S[0] + ext_S_cell[c]);
    qx = q_c_tmp[0][c] + S[1];
    qy = q_c_tmp[1][c] + S[2];

    f_temp0[0][c] = h;
    f_temp1[0][c] = qx;
    f_temp1[1][c] = qy;
  }
}

} // namespace ShallowWater
} // namespace Amanzi
