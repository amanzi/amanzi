/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <set>

#include "errors.hh"
#include "OperatorDefs.hh"

#include "Flow_PK.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* TODO: Verify that a BC has been applied to every boundary face.
* Right now faces without BC are considered no-mass-flux.
****************************************************************** */
void Flow_PK::VV_ValidateBCs() const
{
  // Create sets of the face indices belonging to each BC type.
  std::set<int> pressure_faces, head_faces, flux_faces;

  for (int i =0; i < bcs_.size(); i++) {
    if (bcs_[i]->bc_name() == "pressure") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        pressure_faces.insert(it->first);
      }  
    }

    if (bcs_[i]->bc_name() == "flux") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        flux_faces.insert(it->first);
      }
    }

    if (bcs_[i]->bc_name() == "head") {
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        head_faces.insert(it->first);
      }
    }
  }

  std::set<int> overlap;
  std::set<int>::iterator overlap_end;
  int local_overlap, global_overlap;

  // Check for overlap between pressure and static head BC.
  std::set_intersection(pressure_faces.begin(), pressure_faces.end(),
                        head_faces.begin(), head_faces.end(),
                        std::inserter(overlap, overlap.end()));
  local_overlap = overlap.size();
  mesh_->get_comm()->SumAll(&local_overlap, &global_overlap, 1);  // this will over count ghost faces

  if (global_overlap != 0) {
    Errors::Message msg;
    std::stringstream s;
    s << global_overlap;
    msg << "Flow PK: static head BC overlap Dirichlet BC on "
        << s.str().c_str() << " faces\n";
    Exceptions::amanzi_throw(msg);
  }

  // Check for overlap between pressure and flux BC.
  overlap.clear();
  std::set_intersection(pressure_faces.begin(), pressure_faces.end(),
                        flux_faces.begin(), flux_faces.end(),
                        std::inserter(overlap, overlap.end()));
  local_overlap = overlap.size();
  mesh_->get_comm()->SumAll(&local_overlap, &global_overlap, 1);  // this will over count ghost faces

  if (global_overlap != 0) {
    Errors::Message msg;
    std::stringstream s;
    s << global_overlap;
    msg << "Flow PK: flux BC overlap Dirichlet BC on "
        << s.str().c_str() << " faces\n";
    Exceptions::amanzi_throw(msg);
  }

  // Check for overlap between static head and flux BC.
  overlap.clear();
  std::set_intersection(head_faces.begin(), head_faces.end(),
                        flux_faces.begin(), flux_faces.end(),
                        std::inserter(overlap, overlap.end()));
  local_overlap = overlap.size();
  mesh_->get_comm()->SumAll(&local_overlap, &global_overlap, 1);  // this will over count ghost faces

  if (global_overlap != 0) {
    Errors::Message msg;
    std::stringstream s;
    s << global_overlap;
    msg << "Flow PK: flux BC overlap static head BC on "
        << s.str().c_str() << " faces\n";
    Exceptions::amanzi_throw(msg);
  }
}


/* *******************************************************************
* Reports water balance.
******************************************************************* */
void Flow_PK::VV_ReportWaterBalance(const Teuchos::Ptr<State>& S) const
{
  const Epetra_MultiVector& phi = *S->GetFieldData(porosity_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& flux = *S->GetFieldData(darcy_flux_key_)->ViewComponent("face", true);
  const Epetra_MultiVector& ws = *S->GetFieldData(saturation_liquid_key_)->ViewComponent("cell", false);

  std::vector<int>& bc_model = op_bc_->bc_model();
  double mass_bc_dT = WaterVolumeChangePerSecond(bc_model, flux) * rho_ * dt_;

  double mass_amanzi = 0.0;
  for (int c = 0; c < ncells_owned; c++) {
    mass_amanzi += ws[0][c] * rho_ * phi[0][c] * mesh_->cell_volume(c);
  }

  double mass_amanzi_tmp = mass_amanzi, mass_bc_tmp = mass_bc_dT;
  mesh_->get_comm()->SumAll(&mass_amanzi_tmp, &mass_amanzi, 1);
  mesh_->get_comm()->SumAll(&mass_bc_tmp, &mass_bc_dT, 1);

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
void Flow_PK::VV_ReportSeepageOutflow(const Teuchos::Ptr<State>& S, double dT) const
{
  const Epetra_MultiVector& flux = *S->GetFieldData(darcy_flux_key_)->ViewComponent("face");

  int dir, f, c, nbcs(0);
  double tmp, outflow(0.0);

  for (int i = 0; i < bcs_.size(); ++i) {
    if (bcs_[i]->bc_name() == "seepage") {
      nbcs++;
      for (auto it = bcs_[i]->begin(); it != bcs_[i]->end(); ++it) {
        f = it->first;
        if (f < nfaces_owned) {
          c = BoundaryFaceGetCell(f);
          const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, c, &dir);
          tmp = flux[0][f] * dir;
          if (tmp > 0.0) outflow += tmp;
        }
      }
    }
  }

  tmp = outflow;
  mesh_->get_comm()->SumAll(&tmp, &outflow, 1);

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
void Flow_PK::VV_PrintHeadExtrema(const CompositeVector& pressure) const
{
  std::vector<int>& bc_model = op_bc_->bc_model();
  std::vector<double>& bc_value = op_bc_->bc_value();

  int flag(0);
  double hmin(1.4e+9), hmax(-1.4e+9);  // diameter of the Sun
  double rho_g = rho_ * fabs(gravity_[dim - 1]);

  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_model[f] == Operators::OPERATOR_BC_DIRICHLET) {
      double z = mesh_->face_centroid(f)[dim - 1]; 
      double h = z + (bc_value[f] - atm_pressure_) / rho_g;
      hmax = std::max(hmax, h);
      hmin = std::min(hmin, h);
      flag = 1;
    }
  }
  int flag_tmp(flag);
  mesh_->get_comm()->MinAll(&flag_tmp, &flag, 1);
  
  double tmp;
  Teuchos::OSTab tab = vo_->getOSTab();

  if (flag) {
    tmp = hmin;  // global extrema
    mesh_->get_comm()->MinAll(&tmp, &hmin, 1);
    tmp = hmax;
    mesh_->get_comm()->MaxAll(&tmp, &hmax, 1);

    *vo_->os() << "boundary head (BCs): min=" << units_.OutputLength(hmin) 
               << ", max=" << units_.OutputLength(hmax) << std::endl;
  }

  // process cell-based quantaties
  const Epetra_MultiVector& pcells = *pressure.ViewComponent("cell");
  double vmin(1.4e+9), vmax(-1.4e+9);
  for (int c = 0; c < ncells_owned; c++) {
    double z = mesh_->cell_centroid(c)[dim - 1];              
    double h = z + (pcells[0][c] - atm_pressure_) / rho_g;
    vmax = std::max(vmax, h);
    vmin = std::min(vmin, h);
  }
  tmp = vmin;  // global extrema
  mesh_->get_comm()->MinAll(&tmp, &vmin, 1);
  tmp = vmax;
  mesh_->get_comm()->MaxAll(&tmp, &vmax, 1);
  *vo_->os() << "domain head (cells): min=" << units_.OutputLength(vmin) 
             << ", max=" << units_.OutputLength(vmax) << std::endl;

  // process face-based quantaties (if any)
  if (pressure.HasComponent("face")) {
    const Epetra_MultiVector& pface = *pressure.ViewComponent("face");
    const auto& fmap = *pressure.Map().Map("face", false);

    for (int f = 0; f < nfaces_owned; f++) {
      double z = mesh_->face_centroid(f)[dim - 1];              
      int g = fmap.FirstPointInElement(f);
      double h = z + (pface[0][g] - atm_pressure_) / rho_g;
      vmax = std::max(vmax, h);
      vmin = std::min(vmin, h);
    }
    tmp = vmin;  // global extrema
    mesh_->get_comm()->MinAll(&tmp, &vmin, 1);
    tmp = vmax;
    mesh_->get_comm()->MaxAll(&tmp, &vmax, 1);
    *vo_->os() << "domain head (cells + faces): min=" << units_.OutputLength(vmin) 
               << ", max=" << units_.OutputLength(vmax) << std::endl;
  }
}

 
/* ******************************************************************
* Calculate source extrema.                                   
****************************************************************** */
void Flow_PK::VV_PrintSourceExtrema() const
{
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    double smin(1.0e+99), smax(0.0);
    std::vector<double> volumes;

    for (int i = 0; i < srcs.size(); ++i) {
      for (auto it = srcs[i]->begin(); it != srcs[i]->end(); ++it) {
        smin = std::min(smin, it->second[0]);
        smax = std::max(smax, it->second[0]);
      }
      volumes.push_back(srcs[i]->domain_volume());
    }

    double tmp(smin);
    mesh_->get_comm()->MinAll(&tmp, &smin, 1);
    tmp = smax;
    mesh_->get_comm()->MaxAll(&tmp, &smax, 1);

    if (MyPID == 0) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *vo_->os() << "sources: min=" << smin << " max=" << smax << "  volumes: ";
      for (int i = 0; i < std::min(5, (int)srcs.size()); ++i) {
        *vo_->os() << volumes[i] << " ";
      }
      *vo_->os() << std::endl;
    }
  }
}


/* ****************************************************************
* 
**************************************************************** */
void Flow_PK::OutputTimeHistory(
    const Teuchos::ParameterList& plist, std::vector<dt_tuple>& dT_history)
{
  if (plist.isParameter("plot time history") && 
      vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "saving time history in file flow_dt_history.txt..." << std::endl;

    char file_name[30];
    sprintf(file_name, "flow_dt_history_%d.txt", ti_phase_counter++);

    std::ofstream ofile;
    ofile.open(file_name);

    for (double n = 0; n < dT_history.size(); n++) {
      ofile << std::setprecision(10) << dT_history[n].first / FLOW_YEAR << " " << dT_history[n].second << std::endl;
    }
    ofile.close();
  }
}

}  // namespace Flow
}  // namespace Amanzi

