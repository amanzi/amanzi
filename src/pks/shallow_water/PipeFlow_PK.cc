/*
  Pipe Flow PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Pipe flow model inherited from shallow water model.

  Author: Giacomo Capodaglio (gcapodaglio@lanl.gov)
*/

#include "PipeFlow_PK.hh"
#include "WaterDepthEvaluator.hh"
#include "PressureHeadEvaluator.hh"
#include "NumericalFluxFactory.hh"
#include "PK_DomainFunctionFactory.hh"

namespace Amanzi {
namespace ShallowWater {

using CV_t = CompositeVector;
using CVS_t = CompositeVectorSpace;

PipeFlow_PK::PipeFlow_PK(Teuchos::ParameterList& pk_tree,
                         const Teuchos::RCP<Teuchos::ParameterList>& glist,
                         const Teuchos::RCP<State>& S,
                         const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, glist, S, soln), ShallowWater_PK(pk_tree, glist, S, soln)
{
  // primary variables
  primary_variable_key_ = Keys::getKey(domain_, "wetted_area");
  prev_primary_variable_key_ = Keys::getKey(domain_, "prev_wetted_area");
  wetted_angle_key_ = Keys::getKey(domain_, "wetted_angle");
  water_depth_key_ = Keys::getKey(domain_, "water_depth");
  direction_key_ = sw_list_->get<std::string>("direction key", "");

  // other parameters
  Manning_coeff_ = sw_list_->get<double>("Manning coefficient", 0.005);
  celerity_ = sw_list_->get<double>("celerity", 2); // m/s
  diameter_key_ = Keys::readKey(*sw_list_, domain_, "diameter", "diameter");

  Teuchos::ParameterList vlist;
  vlist.sublist("verbose object") = sw_list_->sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("PipeFlow", vlist));
}

//--------------------------------------------------------------
// Register fields and field evaluators with the state
// Conservative variables: (A, Au)
//--------------------------------------------------------------
void
PipeFlow_PK::Setup()
{
  ShallowWater_PK::Setup();

  pressure_head_key_ = Keys::getKey(domain_, "pressure_head");

  // -- wetted angle
  if (!S_->HasRecord(wetted_angle_key_)) {
    S_->Require<CV_t, CVS_t>(wetted_angle_key_, Tags::DEFAULT, passwd_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    AddDefaultPrimaryEvaluator(S_, wetted_angle_key_);
  }

  // -- pipe diameter
  if (!S_->HasRecord(diameter_key_)) {
    S_->Require<CV_t, CVS_t>(diameter_key_, Tags::DEFAULT, diameter_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    if (!diameter_key_.empty()) {
      S_->RequireEvaluator(diameter_key_, Tags::DEFAULT);
    } else {
      Errors::Message msg;
      msg << "Pipe diameter needs to be specified \n";
      Exceptions::amanzi_throw(msg);
    }
  }

  // -- water depth
  if (!S_->HasRecord(water_depth_key_)) {
    S_->Require<CV_t, CVS_t>(water_depth_key_, Tags::DEFAULT, water_depth_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    Teuchos::ParameterList elist(water_depth_key_);
    elist.set<std::string>("my key", water_depth_key_).set<std::string>("tag", "");
    auto eval = Teuchos::rcp(new WaterDepthEvaluator(elist));
    S_->SetEvaluator(water_depth_key_, Tags::DEFAULT, eval);
  }

  // -- pressure head
  if (!S_->HasRecord(pressure_head_key_)) {
    S_->Require<CV_t, CVS_t>(pressure_head_key_, Tags::DEFAULT, pressure_head_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    Teuchos::ParameterList elist(pressure_head_key_);
    elist.set<std::string>("my key", pressure_head_key_)
      .set<std::string>("tag", Tags::DEFAULT.get())
      .set<double>("celerity", celerity_);
    auto eval = Teuchos::rcp(new PressureHeadEvaluator(elist));
    S_->SetEvaluator(pressure_head_key_, Tags::DEFAULT, eval);
  }

  // -- direction
  if (!S_->HasRecord(direction_key_)) {
    S_->Require<CV_t, CVS_t>(direction_key_, Tags::DEFAULT, direction_key_)
      .SetMesh(mesh_)
      ->SetGhosted(true)
      ->SetComponent("cell", AmanziMesh::CELL, 2);

    if (direction_key_.empty()) {
      Errors::Message msg;
      msg << "Pipe direction needs to be specified \n";
      Exceptions::amanzi_throw(msg);

    } else {
      S_->RequireEvaluator(direction_key_, Tags::DEFAULT);
    }
  }
}

//--------------------------------------------------------------------
// Discretization of the friction source term
//--------------------------------------------------------------------
double
PipeFlow_PK::NumericalSourceFriction(double h,
                                     double qx,
                                     double qy,
                                     double WettedAngle,
                                     int component,
                                     double PipeD)
{
  double S1 = 0.0;
  double num = 0.0;

  //we have to raise this to the power of 7/3 below so the tolerance needs to be stricter
  if (std::fabs(h) > 1.e-10) {
    if (WettedAngle >= 0.0) {
      double WettedPerimeter = 0.5 * PipeD * WettedAngle;
      num = -g_ * Manning_coeff_ * Manning_coeff_ * pow(WettedPerimeter, 4.0 / 3.0) *
            std::fabs(qx) * qx;
      double denom = pow(h, 7.0 / 3.0);
      S1 = num / denom;
    } else { //junction
      if (component == 0) {
        num = -g_ * Manning_coeff_ * Manning_coeff_ * qx * sqrt(qx * qx + qy * qy);
      } else {
        num = -g_ * Manning_coeff_ * Manning_coeff_ * qy * sqrt(qx * qx + qy * qy);
      }
      double denom = pow(h, 7.0 / 3.0);
      S1 = num / denom;
    }
  }

  return S1;
}

//--------------------------------------------------------------------
// Discretization of the bed slope source term
//--------------------------------------------------------------------
std::vector<double>
PipeFlow_PK::NumericalSourceBedSlope(int c,
                                     double htc,
                                     double Bc,
                                     double Bmax,
                                     const Epetra_MultiVector& B_n,
                                     double PipeD,
                                     std::vector<int> bc_model,
                                     std::vector<double> bc_value_h)
{
  std::vector<double> S(3, 0.0);
  std::vector<double> V_rec(2, 0.0);
  std::vector<double> UL(3), UR(3);
  UL[2] = PipeD;
  UR[2] = PipeD;

  if (std::fabs(htc - Bc) >= 1.e-15) { //cell is not dry

    int dir;
    auto cfaces = mesh_->getCellFaces(c);
    double vol = mesh_->getCellVolume(c);

    double BGrad = 0.0;
    double OtherTermLeft = 0.0;
    double OtherTermRight = 0.0;
    double FaceAreaL = 0.0;
    double FaceAreaR = 0.0;
    double denomL = 0.0;
    double denomR = 0.0;
    double tol = 1.e-9;

    for (int n = 0; n < cfaces.size(); ++n) {
      int f = cfaces[n];
      double farea = mesh_->getFaceArea(f);
      AmanziGeometry::Point normal = mesh_->getFaceNormal(f, c, &dir);
      ProjectNormalOntoMeshDirection(c, normal);

      if (normal[0] > (farea * tol)) { //this identifies the j+1/2 face

        auto cells = mesh_->getFaceCells(f);
        int c1 = cells[0];
        int c2 = (cells.size() == 2) ? cells[1] : -1;
        if (c1 > ncells_owned_ && c2 == -1) continue;
        if (c2 > ncells_owned_) std::swap(c1, c2);

        BGrad += BathymetryEdgeValue(f, B_n); //B_(j+1/2)

        V_rec = ComputeFieldsOnEdge(c1, f, htc, Bc, Bmax, B_n, PipeD);

        UL[0] = V_rec[0]; //wetted area
        UL[1] = V_rec[1]; //wetted angle

        FaceAreaL = mesh_->getFaceArea(f);
        denomL = TotalDepthEdgeValue(c1, f, htc, Bc, Bmax, B_n) - BathymetryEdgeValue(f, B_n);

      }

      else if (normal[0] < -(farea * tol)) { //this identifies the j-1/2 face

        BGrad -= BathymetryEdgeValue(f, B_n); //B_(j+1/2) - //B_(j-1/2)

        auto cells = mesh_->getFaceCells(f);
        int c1 = cells[0];
        int c2 = (cells.size() == 2) ? cells[1] : -1;
        if (c1 > ncells_owned_ && c2 == -1) continue;
        if (c2 > ncells_owned_) std::swap(c1, c2);

        if (c2 == -1) {
          if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
            UR[0] = bc_value_h[f];
            UR[1] = ComputeWettedAngleNewton(bc_value_h[f], PipeD);
          } else {
            // default outflow BC
            V_rec = ComputeFieldsOnEdge(c1, f, htc, Bc, Bmax, B_n, PipeD);

            UR[0] = V_rec[0]; //wetted area
            UR[1] = V_rec[1]; //wetted angle
          }
          denomR = TotalDepthEdgeValue(c1, f, htc, Bc, Bmax, B_n) - BathymetryEdgeValue(f, B_n);

        } else {
          V_rec = ComputeFieldsOnEdge(c2, f, htc, Bc, Bmax, B_n, PipeD);

          UR[0] = V_rec[0];
          UR[1] = V_rec[1];

          denomR = TotalDepthEdgeValue(c2, f, htc, Bc, Bmax, B_n) - BathymetryEdgeValue(f, B_n);
        }

        FaceAreaR = mesh_->getFaceArea(f);
      }
    }

    OtherTermLeft = ComputeHydrostaticPressureForce(UL);
    OtherTermRight = ComputeHydrostaticPressureForce(UR);

    BGrad /= vol; //(B_(j+1/2) - B_(j-1/2)) / cellVol

    double denom = denomL - denomR;
    S[1] = std::fabs(BGrad) < 1.e-14 ?
             0.0 :
             -(FaceAreaL * OtherTermLeft - FaceAreaR * OtherTermRight) * BGrad / denom;

  } // closes cell is not dry

  return S;
}

//--------------------------------------------------------------------
// Compute pressure head
//--------------------------------------------------------------------
double
PipeFlow_PK::ComputePressureHead_(double WettedArea, double PipeD)
{
  double PipeCrossSection = M_PI * 0.25 * PipeD * PipeD;
  return (celerity_ * celerity_ * (WettedArea - PipeCrossSection)) / (g_ * PipeCrossSection);
}

//--------------------------------------------------------------------
// Discretization of the bed slope source term for junction
//--------------------------------------------------------------------
std::vector<double>
PipeFlow_PK::NumericalSourceBedSlope(int c,
                                     double htc,
                                     double Bc,
                                     double Bmax,
                                     const Epetra_MultiVector& B_n,
                                     double PipeD)
{
  std::vector<double> S(3, 0.0);
  std::vector<double> V_rec(2, 0.0);
  std::vector<double> UL(3), UR(3);
  UL[2] = PipeD;
  UR[2] = PipeD;

  if (std::fabs(htc - Bc) >= 1.e-15) { //cell is not dry

    auto cfaces = mesh_->getCellFaces(c);
    double vol = mesh_->getCellVolume(c);

    std::vector<int> xMaxFace(cfaces.size(), 0);
    std::vector<int> xMinFace(cfaces.size(), 0);
    std::vector<int> yMaxFace(cfaces.size(), 0);
    std::vector<int> yMinFace(cfaces.size(), 0);

    double tmp = -1.0e+10;
    int indexToSave = 0;
    for (int n = 0; n < cfaces.size(); ++n) {
      const auto& xf = mesh_->getFaceCentroid(cfaces[n]);
      if (xf[0] >= tmp) {
        tmp = xf[0];
        indexToSave = n;
      }
    }
    xMaxFace[indexToSave] = 1;

    tmp = 1.e+10;
    indexToSave = 0;
    for (int n = 0; n < cfaces.size(); ++n) {
      const auto& xf = mesh_->getFaceCentroid(cfaces[n]);
      if (xf[0] <= tmp) {
        tmp = xf[0];
        indexToSave = n;
      }
    }
    xMinFace[indexToSave] = 1;

    tmp = -1.0e+10;
    indexToSave = 0;
    for (int n = 0; n < cfaces.size(); ++n) {
      const auto& xf = mesh_->getFaceCentroid(cfaces[n]);
      if (xf[1] >= tmp) {
        tmp = xf[1];
        indexToSave = n;
      }
    }
    yMaxFace[indexToSave] = 1;

    tmp = 1.e+10;
    indexToSave = 0;
    for (int n = 0; n < cfaces.size(); ++n) {
      const auto& xf = mesh_->getFaceCentroid(cfaces[n]);
      if (xf[1] <= tmp) {
        tmp = xf[1];
        indexToSave = n;
      }
    }
    yMinFace[indexToSave] = 1;

    double BGrad = 0.0;
    double OtherTermLeft = 0.0;
    double OtherTermRight = 0.0;
    double denomL = 0.0;
    double denomR = 0.0;
    int c1, c2;
    for (int n = 0; n < cfaces.size(); ++n) {
      int f = cfaces[n];
      const auto& xf = mesh_->getFaceCentroid(f);
      auto cells = mesh_->getFaceCells(f);
      c1 = cells[0];
      c2 = (cells.size() == 2) ? cells[1] : -1;
      if (c1 > ncells_owned_ && c2 == -1) continue;
      if (c2 > ncells_owned_) std::swap(c1, c2);

      if (xMaxFace[n] == 1 && c2 != -1) { //this identifies the j+1/2 face

        c1 = cells[0];
        c2 = (cells.size() == 2) ? cells[1] : -1;
        if (c1 > ncells_owned_ && c2 == -1) continue;
        if (c2 > ncells_owned_) std::swap(c1, c2);

        BGrad += BathymetryEdgeValue(f, B_n); //B_(j+1/2)

        V_rec = ComputeFieldsOnEdge(c1, f, htc, Bc, Bmax, B_n, PipeD);

        UL[0] = V_rec[0]; //wetted area
        UL[1] = V_rec[1]; //wetted angle

        denomL = TotalDepthEdgeValue(c1, f, htc, Bc, Bmax, B_n) - BathymetryEdgeValue(f, B_n);

      }

      else if (xMinFace[n] == 1 && c2 != -1) { //this identifies the j-1/2 face

        BGrad -= BathymetryEdgeValue(f, B_n); //B_(j+1/2) - //B_(j-1/2)

        V_rec = ComputeFieldsOnEdge(c2, f, htc, Bc, Bmax, B_n, PipeD);

        UR[0] = V_rec[0];
        UR[1] = V_rec[1];

        denomR = TotalDepthEdgeValue(c2, f, htc, Bc, Bmax, B_n) - BathymetryEdgeValue(f, B_n);
      }
    }

    OtherTermLeft = ComputeHydrostaticPressureForce(UL);
    OtherTermRight = ComputeHydrostaticPressureForce(UR);

    BGrad /= vol; //(B_(j+1/2) - B_(j-1/2)) / cellVol

    double denom = denomL - denomR;
    S[1] = std::fabs(BGrad) < 1.e-14 ? 0.0 : -(OtherTermLeft - OtherTermRight) * BGrad / denom;

    BGrad = 0.0;
    OtherTermLeft = 0.0;
    OtherTermRight = 0.0;
    denomL = 0.0;
    denomR = 0.0;

    for (int n = 0; n < cfaces.size(); ++n) {
      int f = cfaces[n];
      const auto& xf = mesh_->getFaceCentroid(f);
      auto cells = mesh_->getFaceCells(f);
      c1 = cells[0];
      c2 = (cells.size() == 2) ? cells[1] : -1;
      if (c1 > ncells_owned_ && c2 == -1) continue;
      if (c2 > ncells_owned_) std::swap(c1, c2);

      if (yMaxFace[n] == 1 && c2 != -1) { //this identifies the j+1/2 face

        BGrad += BathymetryEdgeValue(f, B_n); //B_(j+1/2)

        V_rec = ComputeFieldsOnEdge(c1, f, htc, Bc, Bmax, B_n, PipeD);

        UL[0] = V_rec[0]; //wetted area
        UL[1] = V_rec[1]; //wetted angle

        denomL = TotalDepthEdgeValue(c1, f, htc, Bc, Bmax, B_n) - BathymetryEdgeValue(f, B_n);

      }

      else if (yMinFace[n] == 1 && c2 != -1) { //this identifies the j-1/2 face

        BGrad -= BathymetryEdgeValue(f, B_n); //B_(j+1/2) - //B_(j-1/2)

        V_rec = ComputeFieldsOnEdge(c2, f, htc, Bc, Bmax, B_n, PipeD);

        UR[0] = V_rec[0];
        UR[1] = V_rec[1];

        denomR = TotalDepthEdgeValue(c2, f, htc, Bc, Bmax, B_n) - BathymetryEdgeValue(f, B_n);
      }
    }

    OtherTermLeft = ComputeHydrostaticPressureForce(UL);
    OtherTermRight = ComputeHydrostaticPressureForce(UR);

    BGrad /= vol; //(B_(j+1/2) - B_(j-1/2)) / cellVol

    denom = denomL - denomR;
    S[2] = std::fabs(BGrad) < 1.e-14 ? 0.0 : -(OtherTermLeft - OtherTermRight) * BGrad / denom;

  } // closes cell is not dry

  return S;
}


//--------------------------------------------------------------
// Calculation of timestep limited by the CFL condition
//--------------------------------------------------------------
double
PipeFlow_PK::get_dt()
{
  ComputeCellArrays_();
  int dir;
  double d, d_min = 1.e10, vn, dt = 1.e10, dt_dry = 1.e-1;
  double h, vx, vy;
  double vnMax = 0.0;
  double hMax = 0.0;
  double vnMin = 0.0;
  double hMin = 0.0;

  const auto& h_c = *S_->Get<CV_t>(primary_variable_key_).ViewComponent("cell", true);
  const auto& vel_c = *S_->Get<CV_t>(velocity_key_).ViewComponent("cell", true);

  for (int cell = 0; cell < model_cells_owned_.size(); cell++) {
    int c = model_cells_owned_[cell];
    const Amanzi::AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);

    auto cfaces = mesh_->getCellFaces(c);

    for (int n = 0; n < cfaces.size(); ++n) {
      int f = cfaces[n];
      bool skipFace = false;
      double farea = mesh_->getFaceArea(f);
      const auto& xf = mesh_->getFaceCentroid(f);
      ;
      AmanziGeometry::Point normal = mesh_->getFaceNormal(f, c, &dir);
      normal /= farea;
      ProjectNormalOntoMeshDirection(c, normal);
      SkipFace(normal, skipFace);

      if (!skipFace) {
        h = h_c[0][c];
        vx = vel_c[0][c];
        vy = vel_c[1][c];

        // computing local (cell, face) timestep using Kurganov's estimate d / (2a)
        vn = (vx * normal[0] + vy * normal[1]);
        d = norm(xc - xf);
        d_min = std::min(d_min, d);

        dt = std::min(d / std::max((2.0 * (std::abs(vn) + std::sqrt(g_ * h))), 1.e-12), dt);
        vnMax = std::max(std::abs(vn), vnMax);
        hMax = std::max(h, hMax);
        vnMin = std::min(std::abs(vn), vnMin);
        hMin = std::min(h, hMin);
      }
    }
  }

  for (int cell = 0; cell < junction_cells_owned_.size(); cell++) {
    int c = model_cells_owned_[cell];
    const Amanzi::AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);

    auto cfaces = mesh_->getCellFaces(c);

    for (int n = 0; n < cfaces.size(); ++n) {
      int f = cfaces[n];
      double farea = mesh_->getFaceArea(f);
      const auto& xf = mesh_->getFaceCentroid(f);
      AmanziGeometry::Point normal = mesh_->getFaceNormal(f, c, &dir);
      normal /= farea;

      h = h_c[0][c];
      vx = vel_c[0][c];
      vy = vel_c[1][c];

      // computing local (cell, face) timestep using Kurganov's estimate d / (2a)
      vn = (vx * normal[0] + vy * normal[1]);
      d = norm(xc - xf);
      d_min = std::min(d_min, d);

      dt = std::min(d / std::max((2.0 * (std::abs(vn) + std::sqrt(g_ * h))), 1.e-12), dt);
      vnMax = std::max(std::abs(vn), vnMax);
      hMax = std::max(h, hMax);
      vnMin = std::min(std::abs(vn), vnMin);
      hMin = std::min(h, hMin);
    }
  }

  // reduce dt_min when dt is too large for completely dry conditions (h = 0, qx = 0, qy = 0)
  if (dt >= d_min * 1.e8) { dt = d_min * dt_dry; }

  double dt_min;
  mesh_->getComm()->MinAll(&dt, &dt_min, 1);

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "stable dt = " << dt_min << ", cfl = " << cfl_ << std::endl;
    *vo_->os() << "max abs vel normal to face = " << vnMax << ", max primary variable = " << hMax
               << std::endl;
    *vo_->os() << "min abs vel normal to face = " << vnMin << ", min primary variable = " << hMin
               << std::endl;
  }

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH && iters_ == max_iters_) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "switching from reduced to regular cfl=" << cfl_ << std::endl;
  }

  if (iters_ < max_iters_) {
    return 0.1 * cfl_ * dt_min;
  } else {
    return cfl_ * dt_min;
  }
}


//--------------------------------------------------------------------
// Compute hydrostatic pressure force
//--------------------------------------------------------------------
double
PipeFlow_PK::ComputeHydrostaticPressureForce(std::vector<double> Data)
{
  //Data[0] = wetted area
  //Data[1] = wetted angle
  //Data[2] = pipe diameter

  double I = 0.0;
  double PipeCrossSection = M_PI * 0.25 * Data[2] * Data[2];

  if ((0.0 < Data[0] && Data[0] < PipeCrossSection)) { //flow is ventilated (free-surface)

    I = 3.0 * sin(Data[1] * 0.5) - pow(sin(Data[1] * 0.5), 3) -
        3.0 * (Data[1] * 0.5) * cos(Data[1] * 0.5);
    I = I * g_ * pow(Data[2], 3) / 24.0;

  }

  else if (Data[0] >= PipeCrossSection) { //flow is pressurized

    I = g_ * Data[0] * (ComputePressureHead_(Data[0], Data[2]) + sqrt(Data[0] / M_PI));
  }

  return I;
}

//--------------------------------------------------------------------
// Update wetted angle and total depth for pipe model
//--------------------------------------------------------------------
void
PipeFlow_PK::UpdateSecondaryFields()
{
  auto& PrimaryVar_c =
    *S_->GetW<CV_t>(primary_variable_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
  auto& WettedAngle_c =
    *S_->GetW<CV_t>(wetted_angle_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
  auto& TotalDepth_c =
    *S_->GetW<CV_t>(total_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
  auto& B_c = *S_->GetW<CV_t>(bathymetry_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
  auto& PipeD_c =
    *S_->GetW<CV_t>(diameter_key_, Tags::DEFAULT, diameter_key_).ViewComponent("cell", true);

  for (int c = 0; c < model_cells_owned_.size(); ++c) {
    int cell = model_cells_owned_[c];
    WettedAngle_c[0][cell] = ComputeWettedAngleNewton(PrimaryVar_c[0][cell], PipeD_c[0][cell]);
    TotalDepth_c[0][cell] = ComputeTotalDepth(
      PrimaryVar_c[0][cell], B_c[0][cell], WettedAngle_c[0][cell], PipeD_c[0][cell]);
  }
  for (int c = 0; c < junction_cells_owned_.size(); ++c) {
    int cell = junction_cells_owned_[c];
    WettedAngle_c[0][cell] = -1.0;
    // NOTE: the primary variable for junctions store the water depth
    TotalDepth_c[0][cell] = PrimaryVar_c[0][cell] + B_c[0][cell];
  }

  Teuchos::rcp_dynamic_cast<EvaluatorPrimary<CV_t, CVS_t>>(
    S_->GetEvaluatorPtr(wetted_angle_key_, Tags::DEFAULT))
    ->SetChanged();
}


//--------------------------------------------------------------------
// Compute total depth
//--------------------------------------------------------------------
double
PipeFlow_PK::ComputeTotalDepth(double PrimaryVar,
                               double Bathymetry,
                               double WettedAngle,
                               double PipeD)
{
  double TotalDepth = 0.0;
  double PipeCrossSection = M_PI * 0.25 * PipeD * PipeD;

  if (PrimaryVar >= 0.0 && PrimaryVar < PipeCrossSection) {
    TotalDepth = ComputeWaterDepth(WettedAngle, PipeD) + Bathymetry;
  }

  else if (PrimaryVar >= PipeCrossSection) {
    TotalDepth = PipeD + Bathymetry + ComputePressureHead_(PrimaryVar, PipeD);
  }

  else {
    Errors::Message msg;
    msg << "Wetted area is negative in ComputeTotalDepth \n";
    Exceptions::amanzi_throw(msg);
  }

  return TotalDepth;
}


//--------------------------------------------------------------------
// Compute wetted area and wetted angle at edge location
//--------------------------------------------------------------------
std::vector<double>
PipeFlow_PK::ComputeFieldsOnEdge(int c,
                                 int e,
                                 double htc,
                                 double Bc,
                                 double Bmax,
                                 const Epetra_MultiVector& B_n,
                                 double PipeD)
{
  std::vector<double> V_rec(2, 0.0); //vector to return

  double ht_edge = TotalDepthEdgeValue(c, e, htc, Bc, Bmax, B_n);
  double B_edge = BathymetryEdgeValue(e, B_n);
  double h_edge = std::max((ht_edge - B_edge), 0.0);

  V_rec[1] = ComputeWettedAngle(h_edge, PipeD);
  V_rec[0] = ComputeWettedArea(V_rec[1], PipeD);

  return V_rec;
}


//--------------------------------------------------------------------
// Compute wetted angle given water depth
//--------------------------------------------------------------------
double
PipeFlow_PK::ComputeWettedAngle(double WaterDepth, double PipeD)
{
  if (WaterDepth >= PipeD)
    return 2.0 * M_PI; //if pipe is filled wetted angle is two pi

  else
    return 2.0 * acos(1.0 - 2.0 * WaterDepth / PipeD);
}

//--------------------------------------------------------------------
// Compute wetted area given wetted angle
//--------------------------------------------------------------------
double
PipeFlow_PK::ComputeWettedArea(double WettedAngle, double PipeD)
{
  if (WettedAngle >= (2.0 * M_PI))
    return M_PI * 0.25 * PipeD * PipeD; //if pipe is filled, wetted area is the full cross section

  else
    return PipeD * PipeD * 0.125 * (WettedAngle - sin(WettedAngle));
}

//--------------------------------------------------------------------
// Compute water depth given wetted angle
//--------------------------------------------------------------------
double
PipeFlow_PK::ComputeWaterDepth(double WettedAngle, double PipeD)
{
  if (WettedAngle >= (2.0 * M_PI))
    return PipeD;

  else
    return PipeD * 0.5 * (1.0 - cos(WettedAngle * 0.5));
}

//--------------------------------------------------------------------
// Compute wetted angle given wetted area with Newton's method
//--------------------------------------------------------------------
double
PipeFlow_PK::ComputeWettedAngleNewton(double WettedArea, double PipeD)
{
  double tol = 1.e-15;
  unsigned max_iter = 10000;
  double WettedAngle = 0.0;
  double PipeCrossSection = M_PI * 0.25 * PipeD * PipeD;

  if (std::fabs(WettedArea) < 1.e-15) { //cell is dry
    WettedAngle = 0.0;
  } else if (WettedArea >= PipeCrossSection) { //cell is fully flooded
    WettedAngle = 2.0 * M_PI;
  } else { //cell is partially flooded
    unsigned iter = 0;
    if (std::fabs(WettedAngle) < 1.e-15) WettedAngle = M_PI; // change initial guess if was zero
    double err = WettedAngle - sin(WettedAngle) - 8.0 * WettedArea / (PipeD * PipeD);
    while (iter < max_iter && std::fabs(err) > tol) {
      WettedAngle = WettedAngle - err / (1.0 - cos(WettedAngle));
      err = WettedAngle - sin(WettedAngle) - 8.0 * WettedArea / (PipeD * PipeD);
      iter++;
    }
  }

  return WettedAngle;
}

//--------------------------------------------------------------
// Scatter Master To Ghosted Extra Evaluators
//--------------------------------------------------------------
void
PipeFlow_PK::ScatterMasterToGhostedExtraEvaluators()
{
  S_->Get<CV_t>(wetted_angle_key_).ScatterMasterToGhosted("cell");
  S_->Get<CV_t>(water_depth_key_).ScatterMasterToGhosted("cell");
  S_->Get<CV_t>(pressure_head_key_).ScatterMasterToGhosted("cell");
}

//--------------------------------------------------------------
// Update Extra Evaluators
//--------------------------------------------------------------
void
PipeFlow_PK::UpdateExtraEvaluators()
{
  S_->GetEvaluator(water_depth_key_).Update(*S_, passwd_);
  S_->GetEvaluator(pressure_head_key_).Update(*S_, passwd_);
}

//--------------------------------------------------------------
// Set Primary Variable BC
//--------------------------------------------------------------
void
PipeFlow_PK::SetPrimaryVariableBC(Teuchos::RCP<Teuchos::ParameterList>& bc_list)
{
  // -- wetted area BC
  if (bc_list->isSublist("wetted area")) {
    PK_DomainFunctionFactory<ShallowWaterBoundaryFunction> bc_factory(mesh_, S_);

    Teuchos::ParameterList& tmp_list = bc_list->sublist("wetted area");
    for (auto it = tmp_list.begin(); it != tmp_list.end(); ++it) {
      std::string name = it->first;
      if (tmp_list.isSublist(name)) {
        Teuchos::ParameterList& spec = tmp_list.sublist(name);

        Teuchos::RCP<ShallowWaterBoundaryFunction> bc =
          bc_factory.Create(spec, "wetted area", AmanziMesh::NODE, Teuchos::null);
        bc->set_bc_name("wetted area");
        bc->set_type(WhetStone::DOF_Type::SCALAR);
        bcs_.push_back(bc);
      }
    }
  }
}

//--------------------------------------------------------------
// Initialize Variables
//--------------------------------------------------------------
void
PipeFlow_PK::InitializeFields()
{
  S_->GetEvaluator(diameter_key_).Update(*S_, diameter_key_);

  auto& PrimaryVar_c =
    *S_->GetW<CV_t>(primary_variable_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");
  auto& WettedAngle_c =
    *S_->GetW<CV_t>(wetted_angle_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");
  auto& WaterDepth_c =
    *S_->GetW<CV_t>(water_depth_key_, Tags::DEFAULT, water_depth_key_).ViewComponent("cell");
  auto& PressureHead_c =
    *S_->GetW<CV_t>(pressure_head_key_, Tags::DEFAULT, pressure_head_key_).ViewComponent("cell");
  auto& PipeD_c =
    *S_->GetW<CV_t>(diameter_key_, Tags::DEFAULT, diameter_key_).ViewComponent("cell");
  auto& ht_c = *S_->GetW<CV_t>(total_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");
  auto& B_c = *S_->GetW<CV_t>(bathymetry_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");

  for (int cell = 0; cell < ncells_owned_; cell++) {
    double maxDepth = B_c[0][cell] + PipeD_c[0][cell];
    if (ht_c[0][cell] >= maxDepth) { // cell is pressurized

      double PipeCrossSection = M_PI * 0.25 * PipeD_c[0][cell] * PipeD_c[0][cell];
      PressureHead_c[0][cell] = ht_c[0][cell] - PipeD_c[0][cell] - B_c[0][cell];
      PrimaryVar_c[0][cell] =
        (g_ * PipeCrossSection * PressureHead_c[0][cell]) / (celerity_ * celerity_) +
        PipeCrossSection;
      WettedAngle_c[0][cell] = 2.0 * M_PI;
      WaterDepth_c[0][cell] = PipeD_c[0][cell];
    } else if ((std::fabs(ht_c[0][cell] - B_c[0][cell]) < 1.e-15) ||
               (ht_c[0][cell] < B_c[0][cell])) { //cell is dry

      PrimaryVar_c[0][cell] = 0.0;
      WettedAngle_c[0][cell] = 0.0;
      WaterDepth_c[0][cell] = 0.0;

    }

    else if (ht_c[0][cell] < maxDepth && B_c[0][cell] < ht_c[0][cell]) { //cell is ventilated

      WaterDepth_c[0][cell] = ht_c[0][cell] - B_c[0][cell];
      WettedAngle_c[0][cell] = ComputeWettedAngle(WaterDepth_c[0][cell], PipeD_c[0][cell]);
      PrimaryVar_c[0][cell] = ComputeWettedArea(WettedAngle_c[0][cell], PipeD_c[0][cell]);
    }
  }

  S_->GetRecordW(primary_variable_key_, Tags::DEFAULT, passwd_).set_initialized();
  S_->GetRecordW(wetted_angle_key_, Tags::DEFAULT, passwd_).set_initialized();
  S_->GetRecordW(water_depth_key_, Tags::DEFAULT, water_depth_key_).set_initialized();
  S_->GetRecordW(pressure_head_key_, Tags::DEFAULT, pressure_head_key_).set_initialized();
}

//---------------------------------------------------------------
// Initialize cell array of pipe cells (model cells) and junction
//---------------------------------------------------------------
void
PipeFlow_PK::ComputeCellArrays_()
{
  if (!cellArraysInitDone_) {
    S_->Get<CV_t>(diameter_key_).ScatterMasterToGhosted("cell");
    S_->Get<CV_t>(direction_key_).ScatterMasterToGhosted("cell");
    int ncells_wghost =
      mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);

    auto& PrimaryVar_c =
      *S_->GetW<CV_t>(primary_variable_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
    auto& WettedAngle_c =
      *S_->GetW<CV_t>(wetted_angle_key_, Tags::DEFAULT, passwd_).ViewComponent("cell", true);
    auto& WaterDepth_c =
      *S_->GetW<CV_t>(water_depth_key_, Tags::DEFAULT, water_depth_key_).ViewComponent("cell");
    auto& ht_c = *S_->GetW<CV_t>(total_depth_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");
    auto& B_c = *S_->GetW<CV_t>(bathymetry_key_, Tags::DEFAULT, passwd_).ViewComponent("cell");

    junction_cells_owned_.resize(0);
    model_cells_owned_.resize(0);
    junction_cells_wghost_.resize(0);
    model_cells_wghost_.resize(0);

    for (int c = 0; c < ncells_owned_; c++) {
      // both components of pipe direction zero
      // is code for junction cell
      if (IsJunction_(c)) {
        junction_cells_owned_.push_back(c);
      } else {
        model_cells_owned_.push_back(c);
      }
    }

    for (int c = 0; c < ncells_wghost; c++) {
      // both components of pipe direction zero
      // is code for junction cell
      if (IsJunction_(c)) {
        junction_cells_wghost_.push_back(c);
      } else {
        model_cells_wghost_.push_back(c);
      }
    }

    // here we are overwriting the initialization for
    // junction cells. This should have been done within
    // InitializeFields() but at init we cannot read
    // the dir_c field from file so we need to do it here

    for (int c = 0; c < junction_cells_owned_.size(); c++) {
      int cell = junction_cells_owned_[c];
      PrimaryVar_c[0][cell] = ht_c[0][cell] - B_c[0][cell];
      WettedAngle_c[0][cell] = -1.0;
      WaterDepth_c[0][cell] = ht_c[0][cell] - B_c[0][cell];
    }

    cellArraysInitDone_ = true;
  }
}

//--------------------------------------------------------------
// Checks if a pipe cell is a junction
//--------------------------------------------------------------
bool
PipeFlow_PK::IsJunction_(int c)
{
  bool isJunction(false);

  const auto& dir_c = *S_->Get<CV_t>(direction_key_, Tags::DEFAULT).ViewComponent("cell", true);
  // both components of pipe direction equal to zero is the definition of junction cell
  if (std::fabs(dir_c[0][c]) < 1.e-14 && std::fabs(dir_c[1][c]) < 1.e-14) isJunction = true;

  return isJunction;
}

//--------------------------------------------------------------
// Check if a face needs to be skipped in the flux computation
// (it is skipped if parallel to the flow direction)
// note that this function expects the normal to already
// be rotated according to the pipe direction
//--------------------------------------------------------------
void
PipeFlow_PK::SkipFace(AmanziGeometry::Point normal, bool& skipFace)
{
  if (std::fabs(normal[0]) < 1.e-10) skipFace = true;
}

//--------------------------------------------------------------
// Compute dx
//--------------------------------------------------------------
void
PipeFlow_PK::GetDx(const int& cell, double& dx)
{
  const auto& cfaces = mesh_->getCellFaces(cell);
  AmanziGeometry::Point x1, x2;
  for (int n = 0; n < cfaces.size(); ++n) {
    int f = cfaces[n];
    int dir;
    AmanziGeometry::Point f_normal = mesh_->getFaceNormal(f, cell, &dir);
    ProjectNormalOntoMeshDirection(cell, f_normal);
    if (f_normal[0] > 1.e-08) {
      x1 = mesh_->getFaceCentroid(f);
    } else if (f_normal[0] < -1.e-08) {
      x2 = mesh_->getFaceCentroid(f);
    }
  }
  dx = norm(x1 - x2);
}


//--------------------------------------------------------------
// Compute external forcing on cells
//--------------------------------------------------------------
void
PipeFlow_PK::ComputeExternalForcingOnCells(std::vector<double>& forcing)
{
  for (int i = 0; i < srcs_.size(); ++i) {
    for (auto it = srcs_[i]->begin(); it != srcs_[i]->end(); ++it) {
      int c = it->first;
      double dx;
      GetDx(c, dx);
      forcing[c] = it->second[0] / dx; // [m^2 / s] for pipe
    }
  }
}

//--------------------------------------------------------------
// Project normal onto mesh direction
//--------------------------------------------------------------
void
PipeFlow_PK::ProjectNormalOntoMeshDirection(int c, AmanziGeometry::Point& normal)
{
  // the mesh direction is the pipe direction
  auto& dir_c =
    *S_->GetW<CV_t>(direction_key_, Tags::DEFAULT, direction_key_).ViewComponent("cell", true);
  bool counterClockwise = (dir_c[1][c] < 0.0) ? false : true;
  // first find the angle between pipe direction and x-axis
  // this angle has values between 0 and pi
  double angle = acos(dir_c[0][c] / sqrt(dir_c[0][c] * dir_c[0][c] + dir_c[1][c] * dir_c[1][c]));
  // then rotate the unit vectors of the reference frame according to this angle
  std::vector<double> e1(2);
  std::vector<double> e2(2);
  double e1Tmp1 = 1;
  double e1Tmp2 = 0;
  double e2Tmp1 = 0;
  double e2Tmp2 = 1;
  if (counterClockwise) {
    e1[0] = e1Tmp1 * cos(angle) - e1Tmp2 * sin(angle);
    e1[1] = e1Tmp1 * sin(angle) + e1Tmp2 * cos(angle);
    e2[0] = e2Tmp1 * cos(angle) - e2Tmp2 * sin(angle);
    e2[1] = e2Tmp1 * sin(angle) + e2Tmp2 * cos(angle);
  } else {
    e1[0] = e1Tmp1 * cos(angle) + e1Tmp2 * sin(angle);
    e1[1] = -e1Tmp1 * sin(angle) + e1Tmp2 * cos(angle);
    e2[0] = e2Tmp1 * cos(angle) + e2Tmp2 * sin(angle);
    e2[1] = -e2Tmp1 * sin(angle) + e2Tmp2 * cos(angle);
  }
  double e1Norm = sqrt(e1[0] * e1[0] + e1[1] * e1[1]);
  double e2Norm = sqrt(e2[0] * e2[0] + e2[1] * e2[1]);
  // finally, project the normal on rotated frame
  double nTmp1 = normal[0];
  double nTmp2 = normal[1];
  normal[0] = (nTmp1 * e1[0] + nTmp2 * e1[1]) / e1Norm;
  normal[1] = (nTmp1 * e2[0] + nTmp2 * e2[1]) / e2Norm;
}

} // namespace ShallowWater
} // namespace Amanzi
