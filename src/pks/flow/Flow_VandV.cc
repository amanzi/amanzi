/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

*/

#include <set>

#include "errors.hh"
#include "Mesh_Algorithms.hh"
#include "OperatorDefs.hh"
#include "UniqueLocalIndex.hh"

#include "Flow_PK.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* TBW
****************************************************************** */
void
Flow_PK::VV_FractureConservationLaw() const
{
  if (!coupled_to_matrix_ || fabs(dt_) < 1e+10) return;

  const auto& fracture_flux =
    *S_->Get<CompositeVector>("fracture-volumetric_flow_rate").ViewComponent("face", true);
  const auto& matrix_flux =
    *S_->Get<CompositeVector>("volumetric_flow_rate").ViewComponent("face", true);

  const auto& fracture_map =
    S_->Get<CompositeVector>("fracture-volumetric_flow_rate").Map().Map("face", true);
  const auto& matrix_map = S_->Get<CompositeVector>("volumetric_flow_rate").Map().Map("face", true);

  auto mesh_matrix = S_->GetMesh("domain");

  std::vector<int> dirs;
  double err(0.0), flux_max(0.0);

  for (int c = 0; c < ncells_owned; c++) {
    double flux_sum(0.0);
    auto [faces, fdirs] = mesh_->getCellFacesAndDirections(c);
    for (int i = 0; i < faces.size(); i++) {
      int f = faces[i];
      int g = fracture_map->FirstPointInElement(f);
      int ndofs = fracture_map->ElementSize(f);
      if (ndofs > 1) g += Operators::UniqueIndexFaceToCells(*mesh_, f, c);

      flux_sum += fracture_flux[0][g] * dirs[i];
      flux_max = std::max<double>(flux_max, std::fabs(fracture_flux[0][g]));
    }

    // sum into fluxes from matrix
    auto f = mesh_->getEntityParent(AmanziMesh::Entity_kind::CELL, c);
    auto cells = mesh_matrix->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
    int pos = Operators::UniqueIndexFaceToCells(*mesh_matrix, f, cells[0]);

    for (int j = 0; j != cells.size(); ++j) {
      auto [faces, dirs] = mesh_matrix->getCellFacesAndDirections(cells[j]);

      for (int i = 0; i < faces.size(); i++) {
        if (f == faces[i]) {
          int g = matrix_map->FirstPointInElement(f);
          double fln = matrix_flux[0][g + (pos + j) % 2] * dirs[i];
          flux_sum -= fln;
          flux_max = std::max<double>(flux_max, std::fabs(fln));
          break;
        }
      }
    }

    err = std::max<double>(err, fabs(flux_sum));
  }

  if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "maximum error in conservation law:" << err << " flux_max=" << flux_max
               << std::endl;
  }
}


/* ******************************************************************
* TODO: Verify that a BC has been applied to every boundary face.
* Right now faces without BC are considered no-mass-flux.
****************************************************************** */
void
Flow_PK::VV_ValidateBCs() const
{
  // Create sets of the face indices belonging to each BC type.
  std::set<int> pressure_faces, head_faces, flux_faces;

  for (int i = 0; i < bcs_.size(); i++) {
    if (bcs_[i]->get_bc_name() == "pressure") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        pressure_faces.insert(it->first);
      }
    }

    if (bcs_[i]->get_bc_name() == "flux") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) { flux_faces.insert(it->first); }
    }

    if (bcs_[i]->get_bc_name() == "head") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) { head_faces.insert(it->first); }
    }
  }

  std::set<int> overlap;
  std::set<int>::iterator overlap_end;
  int local_overlap, global_overlap;

  // Check for overlap between pressure and static head BC.
  std::set_intersection(pressure_faces.begin(),
                        pressure_faces.end(),
                        head_faces.begin(),
                        head_faces.end(),
                        std::inserter(overlap, overlap.end()));
  local_overlap = overlap.size();
  mesh_->getComm()->SumAll(&local_overlap, &global_overlap, 1); // this will over count ghost faces

  if (global_overlap != 0) {
    Errors::Message msg;
    std::stringstream s;
    s << global_overlap;
    msg << "Flow PK: static head BC overlap Dirichlet BC on " << s.str().c_str() << " faces\n";
    Exceptions::amanzi_throw(msg);
  }

  // Check for overlap between pressure and flux BC.
  overlap.clear();
  std::set_intersection(pressure_faces.begin(),
                        pressure_faces.end(),
                        flux_faces.begin(),
                        flux_faces.end(),
                        std::inserter(overlap, overlap.end()));
  local_overlap = overlap.size();
  mesh_->getComm()->SumAll(&local_overlap, &global_overlap, 1); // this will over count ghost faces

  if (global_overlap != 0) {
    Errors::Message msg;
    std::stringstream s;
    s << global_overlap;
    msg << "Flow PK: flux BC overlap Dirichlet BC on " << s.str().c_str() << " faces\n";
    Exceptions::amanzi_throw(msg);
  }

  // Check for overlap between static head and flux BC.
  overlap.clear();
  std::set_intersection(head_faces.begin(),
                        head_faces.end(),
                        flux_faces.begin(),
                        flux_faces.end(),
                        std::inserter(overlap, overlap.end()));
  local_overlap = overlap.size();
  mesh_->getComm()->SumAll(&local_overlap, &global_overlap, 1); // this will over count ghost faces

  if (global_overlap != 0) {
    Errors::Message msg;
    std::stringstream s;
    s << global_overlap;
    msg << "Flow PK: flux BC overlap static head BC on " << s.str().c_str() << " faces\n";
    Exceptions::amanzi_throw(msg);
  }
}


/* *******************************************************************
* Reports water balance.
******************************************************************* */
void
Flow_PK::VV_ReportWaterBalance(const Teuchos::Ptr<State>& S) const
{
  const auto& phi = *S->Get<CompositeVector>(porosity_key_).ViewComponent("cell");
  const auto& flowrate = *S->Get<CompositeVector>(vol_flowrate_key_).ViewComponent("face", true);
  const auto& ws = *S->Get<CompositeVector>(saturation_liquid_key_).ViewComponent("cell");

  std::vector<int>& bc_model = op_bc_->bc_model();
  double mass_bc_dT = WaterVolumeChangePerSecond(bc_model, flowrate) * rho_ * dt_;

  double mass_amanzi = 0.0;
  for (int c = 0; c < ncells_owned; c++) {
    mass_amanzi += ws[0][c] * rho_ * phi[0][c] * mesh_->getCellVolume(c);
  }

  double mass_amanzi_tmp = mass_amanzi, mass_bc_tmp = mass_bc_dT;
  mesh_->getComm()->SumAll(&mass_amanzi_tmp, &mass_amanzi, 1);
  mesh_->getComm()->SumAll(&mass_bc_tmp, &mass_bc_dT, 1);

  mass_bc += mass_bc_dT;

  Teuchos::OSTab tab = vo_->getOSTab();
  *vo_->os() << "reservoir water mass=" << units_.OutputMass(mass_amanzi)
             << ", total influx=" << units_.OutputMass(mass_bc) << std::endl;

  // if (report_water_balance) {
  //  if (mass_initial == 0.0) mass_initial = mass_amanzi;
  //  double mass_lost = mass_amanzi - mass_bc - mass_initial;
  //  *vo_->os() << " T=" << S->time() << " water balance error=" << mass_lost << " [kg]" << std::endl;
  // }
}


/* *******************************************************************
* Calculate flow out of the current seepage face.
******************************************************************* */
void
Flow_PK::VV_ReportSeepageOutflow(const Teuchos::Ptr<State>& S, double dT) const
{
  const auto& flowrate = *S->Get<CompositeVector>(vol_flowrate_key_).ViewComponent("face");

  int dir, f, c, nbcs(0);
  double tmp, outflow(0.0);

  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->get_bc_name() == "seepage") {
      nbcs++;
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        f = it->first;
        if (f < nfaces_owned) {
          c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh_, f);
          mesh_->getFaceNormal(f, c, &dir);
          tmp = flowrate[0][f] * dir;
          if (tmp > 0.0) outflow += tmp;
        }
      }
    }
  }

  tmp = outflow;
  mesh_->getComm()->SumAll(&tmp, &outflow, 1);

  outflow *= rho_;
  seepage_mass_ += outflow * dT;

  if (MyPID == 0 && nbcs > 0) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "seepage face: flow=" << outflow << " [kg/s],"
               << " total=" << seepage_mass_ << " [kg]" << std::endl;
  }
}


/* *******************************************************************
* Calculates extrema for hydraulic head.
******************************************************************* */
void
Flow_PK::VV_PrintHeadExtrema(const CompositeVector& pressure) const
{
  std::vector<int>& bc_model = op_bc_->bc_model();
  std::vector<double>& bc_value = op_bc_->bc_value();

  int flag(0);
  double hmin(1.4e+9), hmax(-1.4e+9); // diameter of the Sun
  double rho_g = rho_ * fabs(gravity_[dim - 1]);

  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
      double z = mesh_->getFaceCentroid(f)[dim - 1];
      double h = z + (bc_value[f] - atm_pressure_) / rho_g;
      hmax = std::max(hmax, h);
      hmin = std::min(hmin, h);
      flag = 1;
    }
  }
  int flag_tmp(flag);
  mesh_->getComm()->MinAll(&flag_tmp, &flag, 1);

  double tmp;
  Teuchos::OSTab tab = vo_->getOSTab();

  if (flag == 1) {
    tmp = hmin; // global extrema
    mesh_->getComm()->MinAll(&tmp, &hmin, 1);
    tmp = hmax;
    mesh_->getComm()->MaxAll(&tmp, &hmax, 1);

    *vo_->os() << "boundary head (BCs): min=" << units_.OutputLength(hmin)
               << ", max=" << units_.OutputLength(hmax) << std::endl;
  }

  // process cell-based quantaties
  const Epetra_MultiVector& pcells = *pressure.ViewComponent("cell");
  double vmin(1.4e+9), vmax(-1.4e+9);
  for (int c = 0; c < ncells_owned; c++) {
    double z = mesh_->getCellCentroid(c)[dim - 1];
    double h = z + (pcells[0][c] - atm_pressure_) / rho_g;
    vmax = std::max(vmax, h);
    vmin = std::min(vmin, h);
  }
  tmp = vmin; // global extrema
  mesh_->getComm()->MinAll(&tmp, &vmin, 1);
  tmp = vmax;
  mesh_->getComm()->MaxAll(&tmp, &vmax, 1);
  *vo_->os() << "domain head (cells): min=" << units_.OutputLength(vmin)
             << ", max=" << units_.OutputLength(vmax) << std::endl;

  // process face-based quantaties (if any)
  if (pressure.HasComponent("face")) {
    const Epetra_MultiVector& pface = *pressure.ViewComponent("face");
    const auto& fmap = *pressure.Map().Map("face", false);

    for (int f = 0; f < nfaces_owned; f++) {
      double z = mesh_->getFaceCentroid(f)[dim - 1];
      int g = fmap.FirstPointInElement(f);
      double h = z + (pface[0][g] - atm_pressure_) / rho_g;
      vmax = std::max(vmax, h);
      vmin = std::min(vmin, h);
    }
    tmp = vmin; // global extrema
    mesh_->getComm()->MinAll(&tmp, &vmin, 1);
    tmp = vmax;
    mesh_->getComm()->MaxAll(&tmp, &vmax, 1);
    *vo_->os() << "domain head (cells + faces): min=" << units_.OutputLength(vmin)
               << ", max=" << units_.OutputLength(vmax) << std::endl;
  }
}


/* ******************************************************************
* Calculate source extrema.
****************************************************************** */
void
Flow_PK::VV_PrintSourceExtrema() const
{
  int nsrcs = srcs.size();

  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH && nsrcs > 0) {
    Teuchos::RCP<const Epetra_MultiVector> aperture;
    if (flow_on_manifold_)
      aperture = S_->Get<CompositeVector>(aperture_key_, Tags::DEFAULT).ViewComponent("cell");

    double smin(1.0e+99), smax(-1.0e+99);
    std::vector<double> rates(nsrcs, 0.0), volumes(nsrcs, 0.0), areas(nsrcs, 0.0);

    for (int i = 0; i < nsrcs; ++i) {
      for (auto it = srcs[i]->begin(); it != srcs[i]->end(); ++it) {
        int c = it->first;
        if (c < ncells_owned) {
          double tmp = it->second[0];
          smin = std::min(smin, tmp);
          smax = std::max(smax, tmp);

          double vol = mesh_->getCellVolume(c);

          if (flow_on_manifold_) {
            areas[i] += vol;
            rates[i] += tmp * vol;
            vol *= (*aperture)[0][c];
          } else {
            rates[i] += tmp * vol;
          }
          volumes[i] += vol;
        }
      }
    }

    double tmp1(smin), tmp2(smax);
    std::vector<double> aux1(rates), aux2(volumes);
    mesh_->getComm()->MinAll(&tmp1, &smin, 1);
    mesh_->getComm()->MaxAll(&tmp2, &smax, 1);
    mesh_->getComm()->SumAll(aux1.data(), rates.data(), nsrcs);
    mesh_->getComm()->SumAll(aux2.data(), volumes.data(), nsrcs);
    if (flow_on_manifold_) {
      std::vector<double> aux3(areas);
      mesh_->getComm()->SumAll(aux3.data(), areas.data(), nsrcs);
    }

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "sources: total min/max: " << smin << "/" << smax << std::endl;
    for (int i = 0; i < nsrcs; ++i) {
      if (flow_on_manifold_) {
        *vo_->os() << " src #" << i << ": area=" << areas[i] << " m^2"
                   << ", rate=" << rates[i] << " kg/s"
                   << ", mean aperture=" << volumes[i] / areas[i] << " m" << std::endl;
      } else {
        *vo_->os() << " src #" << i << ": volume=" << volumes[i] << " m^3"
                   << ", rate=" << rates[i] << " kg/s" << std::endl;
      }
    }
  }
}


/* ****************************************************************
*
**************************************************************** */
void
Flow_PK::OutputTimeHistory(const Teuchos::ParameterList& plist, std::vector<dt_tuple>& dT_history)
{
  if (plist.isParameter("plot time history") && vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "saving time history in file flow_dt_history.txt..." << std::endl;

    char file_name[30];
    snprintf(file_name, 30, "flow_dt_history_%d.txt", ti_phase_counter++);

    std::ofstream ofile;
    ofile.open(file_name);

    for (double n = 0; n < dT_history.size(); n++) {
      ofile << std::setprecision(10) << dT_history[n].first / FLOW_YEAR << " "
            << dT_history[n].second << std::endl;
    }
    ofile.close();
  }
}

} // namespace Flow
} // namespace Amanzi
