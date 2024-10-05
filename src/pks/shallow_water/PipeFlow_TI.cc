/*
 Shallow water PK

 Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
 Amanzi is released under the three-clause BSD License.
 The terms of use and "as is" disclaimer for this license are
 provided in the top-level COPYRIGHT file.

 Authors: Giacomo Capodaglio (gcapodaglio@lanl.gov)
          Naren Vohra (vohra@lanl.gov)
          Svetlana Tokareva (tokareva@lanl.gov)
 */

#include "errors.hh"

#include "PipeFlow_PK.hh"

namespace Amanzi {
namespace ShallowWater {

void
PipeFlow_PK::FunctionalTimeDerivative(double t, const TreeVector& A, TreeVector& fun)
{
  const auto& h_temp = *A.SubVector(0)->Data()->ViewComponent("cell", true);
  const auto& q_temp = *A.SubVector(1)->Data()->ViewComponent("cell", true);

  auto& f_temp0 = *fun.SubVector(0)->Data()->ViewComponent("cell");
  auto& f_temp1 = *fun.SubVector(1)->Data()->ViewComponent("cell");

  int ncells_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
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

  auto& WettedAngle_c =
    *S_->GetW<CompositeVector>(wetted_angle_key_, passwd_).ViewComponent("cell", true);
  auto& PipeD_c =
    *S_->GetW<CompositeVector>(diameter_key_, diameter_key_).ViewComponent("cell", true);
  auto& dir_c =
    *S_->GetW<CompositeVector>(direction_key_, direction_key_).ViewComponent("cell", true);

  for (int c = 0; c < model_cells_wghost_.size(); ++c) {
    int cell = model_cells_wghost_[c];
    double factor = inverse_with_tolerance(h_temp[0][cell], cell_area2_max_);
    vel_c[0][cell] = factor * q_temp[0][cell];
    vel_c[1][cell] = factor * q_temp[1][cell];
    ht_c[0][cell] =
      ComputeTotalDepth(h_temp[0][cell], B_c[0][cell], WettedAngle_c[0][cell], PipeD_c[0][cell]);
  }

  for (int c = 0; c < junction_cells_wghost_.size(); ++c) {
    int cell = junction_cells_wghost_[c];
    double factor = inverse_with_tolerance(h_temp[0][cell], cell_area2_max_);
    vel_c[0][cell] = factor * q_temp[0][cell];
    vel_c[1][cell] = factor * q_temp[1][cell];
    ht_c[0][cell] = ShallowWater_PK::ComputeTotalDepth(h_temp[0][cell], B_c[0][cell]);
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

  // NOTE: right now we are assuming that the junction cells cannot be boundary cells
  // so the code below does not take into account a scenario where such a thing happens

  // ponded depth or wetted area BCs
  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->get_bc_name() == "wetted area") { // pipe flow
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
          auto cells = mesh_->getFaceCells(f);
          auto nodes = mesh_->getFaceNodes(f);
          int cell = (cells[0] == -1) ? cells[1] : cells[0];
          int n0 = nodes[0], n1 = nodes[1];
          bc_value_h[f] = (bc_value_hn[n0] + bc_value_hn[n1]) / 2.0;
          bc_value_qx[f] = bc_value_h[f] * it->second[0];
          bc_value_qy[f] = bc_value_h[f] * it->second[1];
          bc_model_scalar[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_value_b[f] = (B_n[0][n0] + B_n[0][n1]) / 2.0;
          double WettedAngle_f = ComputeWettedAngleNewton(bc_value_h[f], PipeD_c[0][cell]);
          bc_value_ht[f] =
            ComputeTotalDepth(bc_value_h[f], bc_value_b[f], WettedAngle_f, PipeD_c[0][cell]);
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
          auto cells = mesh_->getFaceCells(f);
          auto nodes = mesh_->getFaceNodes(f);
          int cell = (cells[0] == -1) ? cells[1] : cells[0];
          int n0 = nodes[0], n1 = nodes[1];
          bc_model_scalar[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_value_h[f] = (bc_value_hn[n0] + bc_value_hn[n1]) / 2.0;
          bc_value_b[f] = (B_n[0][n0] + B_n[0][n1]) / 2.0;
          double WettedAngle_f = ComputeWettedAngleNewton(bc_value_h[f], PipeD_c[0][cell]);
          bc_value_ht[f] =
            ComputeTotalDepth(bc_value_h[f], bc_value_b[f], WettedAngle_f, PipeD_c[0][cell]);
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
          auto cells = mesh_->getFaceCells(f);
          auto nodes = mesh_->getFaceNodes(f);
          int cell = (cells[0] == -1) ? cells[1] : cells[0];
          int n0 = nodes[0], n1 = nodes[1];
          bc_model_scalar[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_value_h[f] = (bc_value_hn[n0] + bc_value_hn[n1]) / 2.0;
          bc_value_b[f] = (B_n[0][n0] + B_n[0][n1]) / 2.0;
          double WettedAngle_f = ComputeWettedAngleNewton(bc_value_h[f], PipeD_c[0][cell]);
          bc_value_ht[f] =
            ComputeTotalDepth(bc_value_h[f], bc_value_b[f], WettedAngle_f, PipeD_c[0][cell]);
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
  // the time step in the last (the second) component of local data vector
  std::vector<double> ext_S_cell(ncells_owned, 0.0);
  ComputeExternalForcingOnCells(ext_S_cell);

  // Shallow water equations have the form
  // U_t + F_x(U) + G_y(U) = S(U)

  int dir, c1, c2, ierr;
  double h, qx, qy, factor, dx;

  std::vector<double> FNum_rot(3, 0.0);       // fluxes
  std::vector<double> FNum_rotTmp(3, 0.0);    // fluxes
  std::vector<double> BedSlopeSource;         // bed slope source
  std::vector<double> FrictionSource(2, 0.0); // friction source
  std::vector<double> UL(3), UR(3);           // local state vectors
  std::vector<double> DL(3), DR(3);           // arrays to compute the hydrostatic pressure forces

  // Simplest flux form
  // U_i^{n+1} = U_i^n - dt/vol * (F_{i+1/2}^n - F_{i-1/2}^n) + dt * S_i

  for (int f = 0; f < nfaces_wghost; ++f) {
    double farea = mesh_->getFaceArea(f);
    const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
    bool c1IsJunction = false;
    bool c2IsJunction = false;
    bool skipFace = false;

    auto cells = mesh_->getFaceCells(f);
    c1 = cells[0];
    c2 = (cells.size() == 2) ? cells[1] : -1;
    if (c1 > ncells_owned && c2 == -1) continue;
    if (c2 > ncells_owned) std::swap(c1, c2);
    int cDiamJnct = c1;

    if (IsJunction(c1)) {
      c1IsJunction = true;
      // this also implies
      // c2 is next to the junction
    }

    if (c2 != -1 && IsJunction(c2)) {
      c2IsJunction = true;
      // this also implies
      // c1 is next to the junction
      // not that junctions are asumed to NOT be boundary cells at the moment
      if (c1IsJunction) {
        if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
          Teuchos::OSTab tab = vo_->getOSTab();
          *vo_->os() << "The junction cell cannot abut other junction cells" << std::endl;
        }
        abort();
      }
    }

    AmanziGeometry::Point normal = mesh_->getFaceNormal(f, c1, &dir);
    AmanziGeometry::Point normalNotRotated = mesh_->getFaceNormal(f, c1, &dir);
    AmanziGeometry::Point normalRotated = mesh_->getFaceNormal(f, c1, &dir);
    normal /= farea;
    normalNotRotated /= farea;
    normalRotated /= farea;

    if (!IsJunction(c1)) {
      // c1 is NOT a junction:
      ProjectNormalOntoMeshDirection(c1, normal);
      ProjectNormalOntoMeshDirection(c1, normalRotated);
      SkipFace(normalRotated, skipFace);
    } else {
      // c1 is a junction:
      // in this case we consider the pipe direction of c2 to be
      // the same as c1 which makes sense since the pipe and
      // junction cells are contiguous
      if (c2 != -1) {
        ProjectNormalOntoMeshDirection(c2, normalRotated);
        SkipFace(normalRotated, skipFace);
        cDiamJnct = c2;
      }
    }

    if (!skipFace) {
      // primary fields at the edge
      // V_rec[0]: primary variable
      // V_rec[1]: wetted angle
      std::vector<double> V_rec(2, 0.0);

      V_rec = ComputeFieldsOnEdge(
        c1, f, ht_c[0][c1], B_c[0][c1], B_max[0][c1], B_n, PipeD_c[0][cDiamJnct]);
      ierr = ErrorDiagnostics_(t, c1, V_rec[0]);
      if (ierr < 0) break;

      double qx_rec = discharge_x_grad_->getValue(c1, xf);
      double qy_rec = discharge_y_grad_->getValue(c1, xf);

      if (c1IsJunction || c2IsJunction) {
        qx_rec = vel_c[0][c1] * V_rec[0];
        qy_rec = vel_c[1][c1] * V_rec[0];
      }

      factor = inverse_with_tolerance(V_rec[0], cell_area2_max_);

      double vx_rec = factor * qx_rec;
      double vy_rec = factor * qy_rec;

      // rotating velocity to the face-based coordinate system
      // note: this implicitly assumes that the velocity and the
      // normal and tangent to the face are expressed with respect
      // to the same reference frame.
      // for pipe, that is the pipe reference frame
      double vn, vt;
      vn = vx_rec * normal[0] + vy_rec * normal[1];
      vt = -vx_rec * normal[1] + vy_rec * normal[0];

      UL[0] = V_rec[0];
      UL[1] = V_rec[0] * vn;
      UL[2] = V_rec[0] * vt;
      DL[0] = V_rec[0];
      DL[1] = V_rec[1];
      DL[2] = PipeD_c[0][cDiamJnct];

      if (c2 == -1) {
        DR[2] = DL[2];
        if (bc_model_scalar[f] == Operators::OPERATOR_BC_DIRICHLET) {
          UR[0] = bc_value_h[f];
          UL[0] = UR[0];
          DR[1] = ComputeWettedAngleNewton(bc_value_h[f], PipeD_c[0][c1]);
          DL[1] = DR[1];
          DR[0] = bc_value_h[f];
          DL[0] = DR[0];
        } else {
          UR[0] = UL[0];
          DR[0] = DL[0];
          DR[1] = DL[1];
        }
        if (bc_model_vector[f] == Operators::OPERATOR_BC_DIRICHLET) {
          if (outward_discharge_flag == true) {
            UR[1] = bc_value_qx
              [f]; // This assumes that the BC value is specified after taking the dot product with the normal
            UR[2] = bc_value_qy[f];
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
        cDiamJnct = (c2IsJunction) ? c1 : c2;
        V_rec = ComputeFieldsOnEdge(
          c2, f, ht_c[0][c2], B_c[0][c2], B_max[0][c2], B_n, PipeD_c[0][cDiamJnct]);
        ierr = ErrorDiagnostics_(t, c2, V_rec[0]);
        if (ierr < 0) break;

        if (!c2IsJunction && !c1IsJunction) {
          qx_rec = discharge_x_grad_->getValue(c2, xf);
          qy_rec = discharge_y_grad_->getValue(c2, xf);

          factor = inverse_with_tolerance(V_rec[0], cell_area2_max_);

          vx_rec = factor * qx_rec;
          vy_rec = factor * qy_rec;

          vn = vx_rec * normal[0] + vy_rec * normal[1];
          vt = -vx_rec * normal[1] + vy_rec * normal[0];
        }

        if (c1IsJunction || c2IsJunction) {
          qx_rec = vel_c[0][c2] * V_rec[0];
          qy_rec = vel_c[1][c2] * V_rec[0];

          factor = inverse_with_tolerance(V_rec[0], cell_area2_max_);

          vx_rec = factor * qx_rec;
          vy_rec = factor * qy_rec;

          if (c1IsJunction) {
            // if c1 is a junction the normal has not been rotated,
            // c2 is a pipe cell so the normal should be rotated
            vn = vx_rec * normalRotated[0] + vy_rec * normalRotated[1];
            vt = -vx_rec * normalRotated[1] + vy_rec * normalRotated[0];
          }

          if (c2IsJunction) {
            // if c2 is a junction, it means c1 is not, hence both normal
            // and normalRotated have been rotated. Hence we need to
            // use normalNotRotated
            vn = vx_rec * normalNotRotated[0] + vy_rec * normalNotRotated[1];
            vt = -vx_rec * normalNotRotated[1] + vy_rec * normalNotRotated[0];
          }
        }

        UR[0] = V_rec[0];
        UR[1] = V_rec[0] * vn;
        UR[2] = V_rec[0] * vt;
        DR[0] = V_rec[0];
        DR[1] = V_rec[1];
        DR[2] = PipeD_c[0][cDiamJnct];
      }

      double HPFL = ComputeHydrostaticPressureForce(DL);
      double HPFR = ComputeHydrostaticPressureForce(DR);

      if (!c1IsJunction && !c2IsJunction) {
        FNum_rot = numerical_flux_->Compute(UL, UR, HPFL, HPFR);
      }
      if (c1IsJunction) {
        FNum_rotTmp = numerical_flux_->Compute(UL, UR, HPFL, HPFR);
        FNum_rot[0] = FNum_rotTmp[0] / farea;
        FNum_rot[1] = FNum_rotTmp[1] / farea;
        FNum_rot[2] = FNum_rotTmp[2] / farea;
      }
      if (c2IsJunction) {
        FNum_rot = numerical_flux_->Compute(UL, UR, HPFL, HPFR);
        FNum_rotTmp[0] = FNum_rot[0] / farea;
        FNum_rotTmp[1] = FNum_rot[1] / farea;
        FNum_rotTmp[2] = FNum_rot[2] / farea;
      }

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
        if (!c1IsJunction && !c2IsJunction) {
          vol = mesh_->getCellVolume(c2);
          factor = farea / vol;
          h_c_tmp[0][c2] += h * factor;
          q_c_tmp[0][c2] += qx * factor;
          q_c_tmp[1][c2] += qy * factor;
        } else {
          h = FNum_rotTmp[0];
          if (c2IsJunction) {
            qx = FNum_rotTmp[1] * normalNotRotated[0] - FNum_rotTmp[2] * normalNotRotated[1];
            qy = FNum_rotTmp[1] * normalNotRotated[1] + FNum_rotTmp[2] * normalNotRotated[0];
          }
          if (c1IsJunction) {
            // if c1 is a junction, c2 is not and so we need to use the rotated normal
            qx = FNum_rotTmp[1] * normalRotated[0] - FNum_rotTmp[2] * normalRotated[1];
            qy = FNum_rotTmp[1] * normalRotated[1] + FNum_rotTmp[2] * normalRotated[0];
          }
          vol = mesh_->getCellVolume(c2);
          factor = farea / vol;
          h_c_tmp[0][c2] += h * factor;
          q_c_tmp[0][c2] += qx * factor;
          q_c_tmp[1][c2] += qy * factor;
        }
      }

    } //if !skipFace
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
  for (int c = 0; c < model_cells_owned_.size(); ++c) {
    int cell = model_cells_owned_[c];

    BedSlopeSource = NumericalSourceBedSlope(cell,
                                             ht_c[0][cell],
                                             B_c[0][cell],
                                             B_max[0][cell],
                                             B_n,
                                             PipeD_c[0][cell],
                                             bc_model_scalar,
                                             bc_value_h);
    FrictionSource[0] = NumericalSourceFriction(h_temp[0][cell],
                                                q_temp[0][cell],
                                                q_temp[1][cell],
                                                WettedAngle_c[0][cell],
                                                0,
                                                PipeD_c[0][cell]);

    h = h_c_tmp[0][cell] + BedSlopeSource[0] + ext_S_cell[cell];
    qx = q_c_tmp[0][cell] + BedSlopeSource[1] + FrictionSource[0];
    qy = 0.0;

    f_temp0[0][cell] = h;
    f_temp1[0][cell] = qx;
    f_temp1[1][cell] = qy;
  }

  // NOTE: we are currently assuming that the manholes cannot be
  // on cells where there is a junction of the pipe network
  // so the term ext_S_cell is currently not computed at junctions cells
  for (int c = 0; c < junction_cells_owned_.size(); ++c) {
    int cell = junction_cells_owned_[c];

    BedSlopeSource = NumericalSourceBedSlope(
      cell, ht_c[0][cell], B_c[0][cell], B_max[0][cell], B_n, PipeD_c[0][cell]);

    FrictionSource[0] = NumericalSourceFriction(h_temp[0][cell],
                                                q_temp[0][cell],
                                                q_temp[1][cell],
                                                WettedAngle_c[0][cell],
                                                0,
                                                PipeD_c[0][cell]);
    FrictionSource[1] = NumericalSourceFriction(h_temp[0][cell],
                                                q_temp[0][cell],
                                                q_temp[1][cell],
                                                WettedAngle_c[0][cell],
                                                1,
                                                PipeD_c[0][cell]);

    h = h_c_tmp[0][cell] + BedSlopeSource[0];
    qx = q_c_tmp[0][cell] + BedSlopeSource[1] + FrictionSource[0];
    qy = q_c_tmp[1][cell] + BedSlopeSource[2] + FrictionSource[1];

    f_temp0[0][cell] = h;
    f_temp1[0][cell] = qx;
    f_temp1[1][cell] = qy;
  }
}

} // namespace ShallowWater
} // namespace Amanzi
