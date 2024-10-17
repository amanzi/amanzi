/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

 Authors: Svetlana Tokareva (tokareva@lanl.gov)
          Giacomo Capodaglio (gcapodaglio@lanl.gov)
          Naren Vohra (vohra@lanl.gov)
*/

/*
 Shallow water PK

 */

#include "errors.hh"

#include "ShallowWater_PK.hh"

namespace Amanzi {
namespace ShallowWater {

void
ShallowWater_PK::FunctionalTimeDerivative(double t, const TreeVector& A, TreeVector& fun)
{
  const auto& h_temp = *A.SubVector(0)->Data()->ViewComponent("cell", true);
  const auto& q_temp = *A.SubVector(1)->Data()->ViewComponent("cell", true);

  auto& f_temp0 = *fun.SubVector(0)->Data()->ViewComponent("cell");
  auto& f_temp1 = *fun.SubVector(1)->Data()->ViewComponent("cell");

  int ncells_wghost =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  int nfaces_wghost =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
  int nnodes_wghost =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::ALL);

  // distribute data to ghost cells
  A.SubVector(0)->Data()->ScatterMasterToGhosted("cell");
  A.SubVector(1)->Data()->ScatterMasterToGhosted("cell");

  // get primary and conservative fields
  const auto& B_c = *S_->Get<CompositeVector>(bathymetry_key_).ViewComponent("cell", true);
  const auto& B_n = *S_->Get<CompositeVector>(bathymetry_key_).ViewComponent("node", true);
  const auto& B_max = *S_->Get<CompositeVector>(bathymetry_key_).ViewComponent("cmax", true);

  auto& ht_c = *S_->GetW<CompositeVector>(total_depth_key_, passwd_).ViewComponent("cell", true);
  auto& vel_c = *S_->GetW<CompositeVector>(velocity_key_, passwd_).ViewComponent("cell", true);
  auto& riemann_f =
    *S_->GetW<CompositeVector>(riemann_flux_key_, passwd_).ViewComponent("face", true);

  for (int c = 0; c < ncells_wghost; ++c) {
    double factor = inverse_with_tolerance(h_temp[0][c], cell_area2_max_);
    vel_c[0][c] = factor * q_temp[0][c];
    vel_c[1][c] = factor * q_temp[1][c];
    ht_c[0][c] = h_temp[0][c] + B_c[0][c];
  }

  // allocate memory for temporary fields
  Epetra_MultiVector h_c_tmp(h_temp);
  Epetra_MultiVector q_c_tmp(q_temp);

  h_c_tmp.PutScalar(0.0);
  q_c_tmp.PutScalar(0.0);

  // update boundary conditions given by [h u v]
  for (int i = 0; i < bcs_.size(); ++i) { bcs_[i]->Compute(t, t); }

  std::vector<int> bc_model_scalar(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<int> bc_model_vector(nfaces_wghost, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value_hn(nnodes_wghost, 0.0);
  std::vector<double> bc_value_h(nfaces_wghost, 0.0);
  std::vector<double> bc_value_b(nfaces_wghost, 0.0);
  std::vector<double> bc_value_ht(nfaces_wghost, 0.0);
  std::vector<double> bc_value_qx(nfaces_wghost, 0.0);
  std::vector<double> bc_value_qy(nfaces_wghost, 0.0);
  bool outward_discharge_flag = false;
  unsigned primary_variable_Dirichlet = 0;

  // ponded depth BC
  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->get_bc_name() == "ponded depth") { // shallow water
      // BC is at nodes
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int n = it->first;
        bc_value_hn[n] = it->second[0];
      }
      primary_variable_Dirichlet = 1;
    }
  }

  // velocity or discharge BCs
  // extract primary variable BC at nodes to ensure well-balancedness
  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->get_bc_name() == "velocity") {
      if (!primary_variable_Dirichlet) {
        if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
          Teuchos::OSTab tab = vo_->getOSTab();
          *vo_->os() << "Dirichlet BCs not set for primary variable" << std::endl;
          *vo_->os() << "Supply Dirichlet BCs for primary variable or" << std::endl;
          *vo_->os() << "alternatively set Dirichlet BCs for discharge only" << std::endl;
        }
        abort();
      } else {
        for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
          int f = it->first;
          bc_model_vector[f] = Operators::OPERATOR_BC_DIRICHLET;
          auto nodes = mesh_->getFaceNodes(f);
          int n0 = nodes[0], n1 = nodes[1];
          bc_model_scalar[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_value_h[f] = (bc_value_hn[n0] + bc_value_hn[n1]) / 2.0;
          bc_value_qx[f] = bc_value_h[f] * it->second[0];
          bc_value_qy[f] = bc_value_h[f] * it->second[1];
          bc_model_scalar[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_value_b[f] = (B_n[0][n0] + B_n[0][n1]) / 2.0;
          bc_value_ht[f] = bc_value_h[f] + bc_value_b[f];
        }
      }
    }
    if (bcs_[i]->get_bc_name() == "discharge") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        bc_model_vector[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value_qx[f] = it->second[0];
        bc_value_qy[f] = it->second[1];
        if (primary_variable_Dirichlet) {
          auto nodes = mesh_->getFaceNodes(f);
          int n0 = nodes[0], n1 = nodes[1];
          bc_model_scalar[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_value_h[f] = (bc_value_hn[n0] + bc_value_hn[n1]) / 2.0;
          bc_value_b[f] = (B_n[0][n0] + B_n[0][n1]) / 2.0;
          bc_value_ht[f] = bc_value_h[f] + bc_value_b[f];
        }
      }
    }

    else if (bcs_[i]->get_bc_name() == "outward discharge") {
      outward_discharge_flag = true;
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        int f = it->first;
        bc_model_vector[f] = Operators::OPERATOR_BC_DIRICHLET;
        bc_value_qx[f] = it->second[0];
        bc_value_qy[f] = it->second[1];
        if (primary_variable_Dirichlet) {
          auto nodes = mesh_->getFaceNodes(f);
          int n0 = nodes[0], n1 = nodes[1];
          bc_model_scalar[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_value_h[f] = (bc_value_hn[n0] + bc_value_hn[n1]) / 2.0;
          bc_value_b[f] = (B_n[0][n0] + B_n[0][n1]) / 2.0;
          bc_value_ht[f] = bc_value_h[f] + bc_value_b[f];
        }
      }
    }
  }

  // limited reconstructions using boundary data
  // total depth
  auto tmp1 =
    S_->GetW<CompositeVector>(total_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
  total_depth_grad_->Compute(tmp1);
  if (use_limiter_)
    limiter_->ApplyLimiter(tmp1, 0, total_depth_grad_, bc_model_scalar, bc_value_ht);
  total_depth_grad_->data()->ScatterMasterToGhosted("cell");

  // additional depth-positivity correction limiting for fully flooded cells
  auto& ht_grad = *total_depth_grad_->data()->ViewComponent("cell", true);

  for (int c = 0; c < ncells_wghost; ++c) {
    auto cnodes = mesh_->getCellNodes(c);

    // calculate maximum bathymetry value on cell nodes
    double Bmax = 0.0;
    for (int i = 0; i < cnodes.size(); ++i) { Bmax = std::max(B_n[0][cnodes[i]], Bmax); }

    // check if cell is fully flooded and proceed with limiting
    if ((ht_c[0][c] > Bmax) && (ht_c[0][c] - B_c[0][c] > 0.0)) {
      auto cfaces = mesh_->getCellFaces(c);

      double alpha = 1.0; // limiter value
      for (int f = 0; f < cfaces.size(); ++f) {
        const auto& xf = mesh_->getFaceCentroid(cfaces[f]);
        double ht_rec = total_depth_grad_->getValue(c, xf);
        if (ht_rec - BathymetryEdgeValue(cfaces[f], B_n) < 0.0) {
          alpha = std::min(alpha,
                           cfl_positivity_ * (BathymetryEdgeValue(cfaces[f], B_n) - ht_c[0][c]) /
                             (ht_rec - ht_c[0][c]));
        }
      }

      ht_grad[0][c] *= alpha;
      ht_grad[1][c] *= alpha;
    }
  }

  // flux
  auto tmp5 = A.SubVector(1)->Data()->ViewComponent("cell", true);
  discharge_x_grad_->Compute(tmp5, 0);
  if (use_limiter_)
    limiter_->ApplyLimiter(tmp5, 0, discharge_x_grad_, bc_model_vector, bc_value_qx);
  discharge_x_grad_->data()->ScatterMasterToGhosted("cell");

  auto tmp6 = A.SubVector(1)->Data()->ViewComponent("cell", true);
  discharge_y_grad_->Compute(tmp6, 1);
  if (use_limiter_)
    limiter_->ApplyLimiter(tmp6, 1, discharge_y_grad_, bc_model_vector, bc_value_qy);
  discharge_y_grad_->data()->ScatterMasterToGhosted("cell");

  // update source (external) terms
  for (int i = 0; i < srcs_.size(); ++i) { srcs_[i]->Compute(t, t); }

  // compute source (external) values
  // coupling submodel="rate" returns volumetric flux [m^3/s] integrated over
  // the timestep in the last (the second) component of local data vector
  std::vector<double> ext_S_cell(ncells_owned_, 0.0);
  ComputeExternalForcingOnCells(ext_S_cell);

  // Shallow water equations have the form
  // U_t + F_x(U) + G_y(U) = S(U)

  int dir, c1, c2, ierr;
  double h, qx, qy, factor;

  std::vector<double> FNum_rot(3, 0.0); // fluxes
  std::vector<double> BedSlopeSource;   // bed slope source
  std::vector<double> UL(3), UR(3);     // local state vectors
  std::vector<double> DL(1), DR(1);     // data to compute the hydrostatic pressure forces

  // Simplest flux form
  // U_i^{n+1} = U_i^n - dt/vol * (F_{i+1/2}^n - F_{i-1/2}^n) + dt * S_i

  for (int f = 0; f < nfaces_wghost; ++f) {
    double farea = mesh_->getFaceArea(f);
    const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);

    auto cells = mesh_->getFaceCells(f);
    c1 = cells[0];
    c2 = (cells.size() == 2) ? cells[1] : -1;
    if (c1 > ncells_owned_ && c2 == -1) continue;
    if (c2 > ncells_owned_) std::swap(c1, c2);

    AmanziGeometry::Point normal = mesh_->getFaceNormal(f, c1, &dir);
    normal /= farea;

    // primary field at the edge
    double V_rec = 0.0;

    V_rec = ComputeFieldOnEdge(c1, f, ht_c[0][c1], B_c[0][c1], B_max[0][c1], B_n);
    ierr = ErrorDiagnostics_(t, c1, V_rec);
    if (ierr < 0) break;

    double qx_rec = discharge_x_grad_->getValue(c1, xf);
    double qy_rec = discharge_y_grad_->getValue(c1, xf);

    factor = inverse_with_tolerance(V_rec, cell_area2_max_);

    double vx_rec = factor * qx_rec;
    double vy_rec = factor * qy_rec;

    // rotating velocity to the face-based coordinate system
    // note: this implicitly assumes that the velocity and the
    // normal and tangent to the face are expressed with respect
    // to the same reference frame.
    // for SW, that is the standard reference frame
    double vn, vt;
    vn = vx_rec * normal[0] + vy_rec * normal[1];
    vt = -vx_rec * normal[1] + vy_rec * normal[0];

    UL[0] = V_rec;
    UL[1] = V_rec * vn;
    UL[2] = V_rec * vt;
    DL[0] = V_rec;

    if (c2 == -1) {
      if (bc_model_scalar[f] == Operators::OPERATOR_BC_DIRICHLET) {
        UR[0] = bc_value_h[f];
        UL[0] = UR[0];
        DR[0] = bc_value_h[f];
        DL[0] = DR[0];
      } else {
        UR[0] = UL[0];
        DR[0] = DL[0];
      }
      if (bc_model_vector[f] == Operators::OPERATOR_BC_DIRICHLET) {
        if (outward_discharge_flag == true) {
          UR[1] = bc_value_qx[f]; // This assumes that BC value is specified by taking the dot product with face normal
          UR[2] = bc_value_qy[f]; // This should probably be 0.
        } else {
          UR[1] = bc_value_qx[f] * normal[0] + bc_value_qy[f] * normal[1];
          UR[2] = -bc_value_qx[f] * normal[1] + bc_value_qy[f] * normal[0];
        }

        UL[1] = UR[1];
        UL[2] = UR[2];
      } else {
        // default outflow BC
        UR[1] = UL[1];
        UR[2] = UL[2];
      }
    } else {
      V_rec = ComputeFieldOnEdge(c2, f, ht_c[0][c2], B_c[0][c2], B_max[0][c2], B_n);
      ierr = ErrorDiagnostics_(t, c2, V_rec);
      if (ierr < 0) break;

      qx_rec = discharge_x_grad_->getValue(c2, xf);
      qy_rec = discharge_y_grad_->getValue(c2, xf);

      factor = inverse_with_tolerance(V_rec, cell_area2_max_);

      vx_rec = factor * qx_rec;
      vy_rec = factor * qy_rec;

      vn = vx_rec * normal[0] + vy_rec * normal[1];
      vt = -vx_rec * normal[1] + vy_rec * normal[0];

      UR[0] = V_rec;
      UR[1] = V_rec * vn;
      UR[2] = V_rec * vt;
      DR[0] = V_rec;
    }

    double HPFL = ComputeHydrostaticPressureForce(DL);
    double HPFR = ComputeHydrostaticPressureForce(DR);

    FNum_rot = numerical_flux_->Compute(UL, UR, HPFL, HPFR);

    h = FNum_rot[0];
    qx = FNum_rot[1] * normal[0] - FNum_rot[2] * normal[1];
    qy = FNum_rot[1] * normal[1] + FNum_rot[2] * normal[0];

    // save riemann mass flux
    riemann_f[0][f] = FNum_rot[0] * farea * dir;

    // add fluxes to temporary fields
    double vol = mesh_->getCellVolume(c1);
    factor = farea / vol;
    h_c_tmp[0][c1] -= h * factor;
    q_c_tmp[0][c1] -= qx * factor;
    q_c_tmp[1][c1] -= qy * factor;

    if (c2 != -1) {
      vol = mesh_->getCellVolume(c2);
      factor = farea / vol;
      h_c_tmp[0][c2] += h * factor;
      q_c_tmp[0][c2] += qx * factor;
      q_c_tmp[1][c2] += qy * factor;
    }
  }

  // sync errors across processors, otherwise any MPI call
  // will lead to an error.
  int ierr_tmp(ierr);
  mesh_->getComm()->MinAll(&ierr_tmp, &ierr, 1);
  if (ierr < 0) {
    Errors::Message msg;
    msg << "Negative primary variable.\n";
    Exceptions::amanzi_throw(msg);
  }

  // sources (bathymetry, flux exchange, etc)
  for (int c = 0; c < ncells_owned_; ++c) {
    BedSlopeSource = NumericalSourceBedSlope(c, ht_c[0][c], B_c[0][c], B_max[0][c], B_n);

    h = h_c_tmp[0][c] + BedSlopeSource[0] + ext_S_cell[c];
    qx = q_c_tmp[0][c] + BedSlopeSource[1];
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
ShallowWater_PK::ErrorDiagnostics_(double t, int c, double h)
{
  if (h < 0.0) {
    const auto& xc = mesh_->getCellCentroid(c);
    if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "negative primary variable in cell " << c << ", xc=(" << xc[0] << ", " << xc[1]
                 << "), primary variable=" << h << std::endl;
    }
    return -1;
  }
  return 0;
}

} // namespace ShallowWater
} // namespace Amanzi
