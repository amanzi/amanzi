/*
 Shallow water PK

 Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
 Amanzi is released under the three-clause BSD License.
 The terms of use and "as is" disclaimer for this license are
 provided in the top-level COPYRIGHT file.

 Author: Svetlana Tokareva (tokareva@lanl.gov)
 */

#include "errors.hh"

#include "ShallowWater_PK.hh"

namespace Amanzi {
namespace ShallowWater {

void
ShallowWater_PK::FunctionalTimeDerivative(double t, const TreeVector& A,
                                          TreeVector& fun)
{

  const auto& h_temp = *A.SubVector(0)->Data()->ViewComponent("cell", true);
  const auto& q_temp = *A.SubVector(1)->Data()->ViewComponent("cell", true);

  auto& f_temp0 = *fun.SubVector(0)->Data()->ViewComponent("cell");
  auto& f_temp1 = *fun.SubVector(1)->Data()->ViewComponent("cell");

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  int ncells_wghost = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
  int nfaces_wghost = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::ALL);
  int nnodes_wghost = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::ALL);

  // distribute data to ghost cells
  A.SubVector(0)->Data()->ScatterMasterToGhosted("cell");
  A.SubVector(1)->Data()->ScatterMasterToGhosted("cell");

  // get primary and conservative fields
  const auto& B_c = *S_->Get<CompositeVector>(bathymetry_key_).ViewComponent("cell", true);
  const auto& B_n = *S_->Get<CompositeVector>(bathymetry_key_).ViewComponent("node", true);
  const auto& B_max = *S_->Get<CompositeVector>(bathymetry_key_).ViewComponent("cmax", true);

  auto& ht_c = *S_->GetW<CompositeVector>(total_depth_key_, passwd_).ViewComponent("cell", true);
  auto& vel_c = *S_->GetW<CompositeVector>(velocity_key_, passwd_).ViewComponent("cell", true);
  auto& riemann_f = *S_->GetW<CompositeVector>(riemann_flux_key_, passwd_).ViewComponent("face", true);

  auto& WettedAngle_c = *S_->GetW<CompositeVector>(wetted_angle_key_, passwd_).ViewComponent("cell", true); 

  for (int c = 0; c < ncells_wghost; ++c) {
    double factor = inverse_with_tolerance(h_temp[0][c], cell_area2_max_);
    vel_c[0][c] = factor * q_temp[0][c];
    vel_c[1][c] = factor * q_temp[1][c];
    if(!hydrostatic_pressure_force_type_){
       ht_c[0][c] = h_temp[0][c] + B_c[0][c];
    }
    else{
       ht_c[0][c] = ComputeTotalDepth(h_temp[0][c], WettedAngle_c[0][c], B_c[0][c]);
    }
  }

  // allocate memory for temporary fields
  Epetra_MultiVector h_c_tmp(h_temp);
  Epetra_MultiVector q_c_tmp(q_temp);

  h_c_tmp.PutScalar(0.0);
  q_c_tmp.PutScalar(0.0);

  // update boundary conditions given by [h u v]
  for (int i = 0; i < bcs_.size(); ++i) { bcs_[i]->Compute(t, t); }

  std::vector<int> bc_model(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value_hn(nnodes_wghost, 0.0);
  std::vector<double> bc_value_h(nfaces_wghost, 0.0);
  std::vector<double> bc_value_b(nfaces_wghost, 0.0);
  std::vector<double> bc_value_ht(nfaces_wghost, 0.0);
  std::vector<double> bc_value_qx(nfaces_wghost, 0.0);
  std::vector<double> bc_value_qy(nfaces_wghost, 0.0);

  // extract velocity and compute qx, qy, h BC at faces
  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->get_bc_name() == "ponded depth") {
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
        int n0 = nodes[0], n1 = nodes[1];

        bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value_h[f] = (bc_value_hn[n0] + bc_value_hn[n1]) / 2.0;
        bc_value_qx[f] = bc_value_h[f] * it->second[0];
        bc_value_qy[f] = bc_value_h[f] * it->second[1];
        bc_value_b[f] = (B_n[0][n0] + B_n[0][n1]) / 2.0;
        if(!hydrostatic_pressure_force_type_){
           bc_value_ht[f] = bc_value_h[f] + bc_value_b[f];
        }
        else{
           double WettedAngle_f = ComputeWettedAngleNewton(bc_value_h[f]);
           bc_value_ht[f] = ComputeTotalDepth(bc_value_h[f], WettedAngle_f, bc_value_b[f]);
        }
      }
    }
  }

  // limited reconstructions using boundary data
  // total depth
  auto tmp1 = S_->GetW<CompositeVector>(total_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
  total_depth_grad_->Compute(tmp1);
  if (use_limiter_)
    limiter_->ApplyLimiter(tmp1, 0, total_depth_grad_, bc_model, bc_value_ht);
  total_depth_grad_->data()->ScatterMasterToGhosted("cell");

  // additional depth-positivity correction limiting for fully flooded cells
  auto& ht_grad = *total_depth_grad_->data()->ViewComponent("cell", true);
  
  for (int c = 0; c < ncells_wghost; ++c ) {
     Amanzi::AmanziMesh::Entity_ID_List cnodes, cfaces;
     mesh_->cell_get_nodes(c, &cnodes);
    
    // calculate maximum bathymetry value on cell nodes
    double Bmax = 0.0;
    for (int i = 0; i < cnodes.size(); ++i) {
      Bmax = std::max(B_n[0][cnodes[i]], Bmax);
    }
    
    // check if cell is fully flooded and proceed with limiting
    if ((ht_c[0][c] > Bmax) && (ht_c[0][c] - B_c[0][c] > 0.0)) {
      Amanzi::AmanziMesh::Entity_ID_List cfaces;
      mesh_->cell_get_faces(c, &cfaces);

      double alpha = 1.0; // limiter value
      for (int f = 0; f < cfaces.size(); ++f) {
        const auto& xf = mesh_->face_centroid(cfaces[f]);
        double ht_rec = total_depth_grad_->getValue(c, xf);
        if (ht_rec - BathymetryEdgeValue(cfaces[f], B_n) < 0.0) {
          alpha = std::min(alpha, cfl_positivity_ * (BathymetryEdgeValue(cfaces[f], B_n) - ht_c[0][c]) / (ht_rec - ht_c[0][c]));
        }
      }
      
      ht_grad[0][c] *= alpha;
      ht_grad[1][c] *= alpha;
    }
  }

  // compute bathymetry gradient for bed slope source
  auto tmp7 = S_->GetW<CompositeVector>(bathymetry_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
  bathymetry_grad_->Compute(tmp7);
  bathymetry_grad_->data()->ScatterMasterToGhosted("cell");

  // flux
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

  int dir, c1, c2, ierr;
  double h, qx, qy, factor;
  AmanziMesh::Entity_ID_List cells;

  std::vector<double> FNum_rot;        // fluxes
  std::vector<double> BedSlopeSource;  // bed slope source
  double FrictionSource = 0.0;         // friction source
  std::vector<double> UL(4), UR(4), U; // local state vectors and wetted angle

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

    std::vector<double> W_rec(2,0.0);
    double ht_rec;
    double B_rec = BathymetryEdgeValue(f, B_n);;
    double h_rec;

    if (!hydrostatic_pressure_force_type_)  {
        ht_rec = TotalDepthEdgeValue(c1, f, ht_c[0][c1], B_c[0][c1], B_max[0][c1], B_n);
        h_rec = ht_rec - B_rec;
        ierr = ErrorDiagnostics_(t, c1, h_rec, B_rec, ht_rec);
        if (ierr < 0) break;
    }
    else{

       W_rec = ComputeWettedQuantitiesEdge(c1, f, ht_c[0][c1], B_c[0][c1], B_max[0][c1], B_n);
       h_rec = W_rec[0];

    }

    double qx_rec = discharge_x_grad_->getValue(c1, xf);
    double qy_rec = discharge_y_grad_->getValue(c1, xf);

    factor = inverse_with_tolerance(h_rec, cell_area2_max_);
  
    double vx_rec = factor * qx_rec;
    double vy_rec = factor * qy_rec;

    // rotating velocity to the face-based coordinate system
    double vn, vt;
    vn = vx_rec * normal[0] + vy_rec * normal[1];
    vt = -vx_rec * normal[1] + vy_rec * normal[0];

    UL[0] = h_rec;
    UL[1] = h_rec * vn;
    UL[2] = h_rec * vt;
    UL[3] = W_rec[1]; 

    if (c2 == -1) {
      if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
        UR[0] = bc_value_h[f];
        UR[1] = bc_value_qx[f] * normal[0] + bc_value_qy[f] * normal[1];
        UR[2] = -bc_value_qx[f] * normal[1] + bc_value_qy[f] * normal[0];
        UR[3] = ComputeWettedAngleNewton(bc_value_h[f]);
      } else {
        // default outflow BC
        UR = UL;
      }
    } else {

      if (!hydrostatic_pressure_force_type_)  {
         ht_rec = TotalDepthEdgeValue(c2, f, ht_c[0][c2], B_c[0][c2], B_max[0][c2], B_n);
         h_rec = ht_rec - B_rec;
         ierr = ErrorDiagnostics_(t, c2, h_rec, B_rec, ht_rec);
         if (ierr < 0) break;
      }
      else{

       W_rec = ComputeWettedQuantitiesEdge(c2, f, ht_c[0][c2], B_c[0][c2], B_max[0][c2], B_n);
       h_rec = W_rec[0];

      }

      qx_rec = discharge_x_grad_->getValue(c2, xf);
      qy_rec = discharge_y_grad_->getValue(c2, xf);
      
      factor = inverse_with_tolerance(h_rec, cell_area2_max_);

      vx_rec = factor * qx_rec;
      vy_rec = factor * qy_rec;

      vn = vx_rec * normal[0] + vy_rec * normal[1];
      vt = -vx_rec * normal[1] + vy_rec * normal[0];

      UR[0] = h_rec;
      UR[1] = h_rec * vn;
      UR[2] = h_rec * vt;
      UR[3] = W_rec[1];
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

  // sync errors across processors, otherwise any MPI call
  // will lead to an error.
  int ierr_tmp(ierr);
  mesh_->get_comm()->MinAll(&ierr_tmp, &ierr, 1);
  if (ierr < 0) {
    Errors::Message msg;
    msg << "Negative ponded depth.\n";
    Exceptions::amanzi_throw(msg);
  }

  // sources (bathymetry, flux exchange, etc)
  // the code should not fail after that
  U.resize(3);
  double ExtraSource;

  for (int c = 0; c < ncells_owned; ++c) {
    U[0] = h_temp[0][c];
    U[1] = q_temp[0][c];
    U[2] = q_temp[1][c];
    U[3] = WettedAngle_c[0][c];

    if (!hydrostatic_pressure_force_type_){
       BedSlopeSource = NumericalSourceBedSlope(c, U[0] + B_c[0][c], B_c[0][c], B_max[0][c], B_n);
       ExtraSource = 1.0;
    }
    else{
       BedSlopeSource = NumericalSourceBedSlope(c, U[0] + B_c[0][c], B_c[0][c], B_max[0][c], B_n, bc_model, bc_value_h);
       ExtraSource = 0.0;
    }
    FrictionSource = NumericalSourceFriction(U[0], U[1], U[3]); 

    h = h_c_tmp[0][c] + (BedSlopeSource[0] + ext_S_cell[c] * ExtraSource);
    qx = q_c_tmp[0][c] + BedSlopeSource[1] + FrictionSource;
    qy = q_c_tmp[1][c] + BedSlopeSource[2];

    f_temp0[0][c] = h;
    f_temp1[0][c] = qx;
    f_temp1[1][c] = qy;
  }

}


//--------------------------------------------------------------
// Error diagnostics
//--------------------------------------------------------------
int
ShallowWater_PK::ErrorDiagnostics_(double t, int c, double h, double B, double ht)
{
  if (h < 0.0) {
    const auto& xc = mesh_->cell_centroid(c);

    if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "negative height in cell " << c << ", xc=(" << xc[0] << ", " << xc[1] << ")"
                 << ", total=" << ht << ", bathymetry=" << B << ", height=" << h << std::endl;
    }
    return -1;
  }
  return 0;
}

} // namespace ShallowWater
} // namespace Amanzi
