/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <set>

#include "errors.hh"
#include "Flow_PK.hh"

namespace Amanzi {
namespace Flow {

/* ****************************************************************
* Construct default state for unit tests. It completes
* initialization of all missed objects in the state.
**************************************************************** */
void Flow_PK::InitializeFields()
{
  // set popular default values
  if (!S_->GetField("porosity", passwd_)->initialized()) {
    S_->GetFieldData("porosity", passwd_)->PutScalar(0.2);
    S_->GetField("porosity", passwd_)->set_initialized();
  }

  if (!S_->GetField("fluid_density", passwd_)->initialized()) {
    *(S_->GetScalarData("fluid_density", passwd_)) = 1000.0;
    S_->GetField("fluid_density", passwd_)->set_initialized();
  }

  if (!S_->GetField("fluid_viscosity", passwd_)->initialized()) {
    *(S_->GetScalarData("fluid_viscosity", passwd_)) = 0.001;
    S_->GetField("fluid_viscosity", passwd_)->set_initialized();
  }

  if (!S_->GetField("gravity", passwd_)->initialized()) {
    Epetra_Vector& gvec = *S_->GetConstantVectorData("gravity", passwd_);
    gvec.PutScalar(0.0);
    gvec[dim - 1] = -9.80;
    S_->GetField("gravity", passwd_)->set_initialized();
  }

  if (!S_->GetField("permeability", passwd_)->initialized()) {
    S_->GetFieldData("permeability", passwd_)->PutScalar(1.0);
    S_->GetField("permeability", passwd_)->set_initialized();
  }

  if (S_->HasField("water_saturation")) {
    if (!S_->GetField("water_saturation", passwd_)->initialized()) {
      S_->GetFieldData("water_saturation", passwd_)->PutScalar(1.0);
      S_->GetField("water_saturation", passwd_)->set_initialized();
    }
  }

  if (S_->HasField("prev_water_saturation")) {
    if (!S_->GetField("prev_water_saturation", passwd_)->initialized()) {
      S_->GetFieldData("prev_water_saturation", passwd_)->PutScalar(1.0);
      S_->GetField("prev_water_saturation", passwd_)->set_initialized();
    }
  }

  if (S_->HasField("specific_storage")) {
    if (!S_->GetField("specific_storage", passwd_)->initialized()) {
      S_->GetFieldData("specific_storage", passwd_)->PutScalar(0.0);
      S_->GetField("specific_storage", passwd_)->set_initialized();
    }
  }

  if (S_->HasField("specific_yield")) {
    if (!S_->GetField("specific_yield", passwd_)->initialized()) {
      S_->GetFieldData("specific_yield", passwd_)->PutScalar(0.0);
      S_->GetField("specific_yield", passwd_)->set_initialized();
    }
  }

  if (!S_->GetField("pressure", passwd_)->initialized()) {
    S_->GetFieldData("pressure", passwd_)->PutScalar(0.0);
    S_->GetField("pressure", passwd_)->set_initialized();
  }

  if (!S_->GetField("hydraulic_head", passwd_)->initialized()) {
    S_->GetFieldData("hydraulic_head", passwd_)->PutScalar(0.0);
    S_->GetField("hydraulic_head", passwd_)->set_initialized();
  }

  if (!S_->GetField("darcy_flux", passwd_)->initialized()) {
    S_->GetFieldData("darcy_flux", passwd_)->PutScalar(0.0);
    S_->GetField("darcy_flux", passwd_)->set_initialized();
  }

  // if (!S_->GetField("darcy_velocity", passwd_)->initialized()) {
  //   S_->GetFieldData("darcy_velocity", passwd_)->PutScalar(0.0);
  //   S_->GetField("darcy_velocity", passwd_)->set_initialized();
  // }
}


/* ******************************************************************
* TODO: Verify that a BC has been applied to every boundary face.
* Right now faces without BC are considered no-mass-flux.
****************************************************************** */
void Flow_PK::VV_ValidateBCs() const
{
  // Create sets of the face indices belonging to each BC type.
  std::set<int> pressure_faces, head_faces, flux_faces;
  Amanzi::Functions::FlowBoundaryFunction::Iterator bc;
  for (bc = bc_pressure->begin(); bc != bc_pressure->end(); ++bc) pressure_faces.insert(bc->first);
  for (bc = bc_head->begin(); bc != bc_head->end(); ++bc) head_faces.insert(bc->first);
  for (bc = bc_flux->begin(); bc != bc_flux->end(); ++bc) flux_faces.insert(bc->first);

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
* Calculates best least square fit for data (h[i], error[i]).                       
******************************************************************* */
void Flow_PK::VV_PrintHeadExtrema(const CompositeVector& pressure) const
{
  double hmin(1.4e+9), hmax(-1.4e+9);  // diameter of the Sun
  double rho_g = rho_ * fabs(gravity_[dim - 1]);
  for (int f = 0; f < nfaces_owned; f++) {
    if (bc_model[f] == FLOW_BC_FACE_PRESSURE) {
      double z = mesh_->face_centroid(f)[dim - 1]; 
      double h = z + (bc_value[f] - atm_pressure_) / rho_g;
      hmax = std::max(hmax, h);
      hmin = std::min(hmin, h);
    }
  }
  double tmp = hmin;  // global extrema
  mesh_->get_comm()->MinAll(&tmp, &hmin, 1);
  tmp = hmax;
  mesh_->get_comm()->MaxAll(&tmp, &hmax, 1);

  Teuchos::OSTab tab = vo_->getOSTab();
  *vo_->os() << "boundary head (BCs): min=" << hmin << " max=" << hmax << " [m]" << std::endl;

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
  *vo_->os() << "domain head (cells): min=" << vmin << " max=" << vmax << " [m]" << std::endl;

  // process face-based quantaties (if any)
  if (pressure.HasComponent("face")) {
    const Epetra_MultiVector& pface = *pressure.ViewComponent("face");

    for (int f = 0; f < nfaces_owned; f++) {
      double z = mesh_->face_centroid(f)[dim - 1];              
      double h = z + (pface[0][f] - atm_pressure_) / rho_g;
      vmax = std::max(vmax, h);
      vmin = std::min(vmin, h);
    }
    tmp = vmin;  // global extrema
    mesh_->get_comm()->MinAll(&tmp, &vmin, 1);
    tmp = vmax;
    mesh_->get_comm()->MaxAll(&tmp, &vmax, 1);
    *vo_->os() << "domain head (cells + faces): min=" << vmin << " max=" << vmax << " [m]" << std::endl;
  }
}

 
/* *******************************************************************
* Calculates best least square fit for data (h[i], error[i]).                       
******************************************************************* */
double bestLSfit(const std::vector<double>& h, const std::vector<double>& error)
{
  double a = 0.0, b = 0.0, c = 0.0, d = 0.0, tmp1, tmp2;

  int n = h.size();
  for (int i = 0; i < n; i++) {
    tmp1 = log(h[i]);
    tmp2 = log(error[i]);
    a += tmp1;
    b += tmp2;
    c += tmp1 * tmp1;
    d += tmp1 * tmp2;
  }

  return (a * b - n * d) / (a * a - n * c);
}


}  // namespace Flow
}  // namespace Amanzi

