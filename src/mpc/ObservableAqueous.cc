/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
      Konstantin Lipnikov
      Daniil Svyatsky
*/

/*
  Multi-Process Coordinator

*/

#include <map>

#include "Evaluator.hh"
#include "RegionPlane.hh"
#include "RegionPolygon.hh"
#include "ReconstructionCellLinear.hh"
#include "Units.hh"

#include "ObservableAqueous.hh"

namespace Amanzi {

/* ******************************************************************
* Delegating constructor
****************************************************************** */
ObservableAqueous::ObservableAqueous(std::string variable,
                                     std::string region,
                                     std::string functional,
                                     Teuchos::ParameterList& plist,
                                     Teuchos::ParameterList& units_plist,
                                     Teuchos::RCP<const AmanziMesh::Mesh> mesh)
  : Observable(variable, region, functional, plist, units_plist, mesh){};


/* ******************************************************************
* Defines
****************************************************************** */
int
ObservableAqueous::ComputeRegionSize()
{
  //int mesh_block_size;

  // check if observation is planar
  obs_planar_ = false;

  Teuchos::RCP<const AmanziGeometry::GeometricModel> gm_ptr = mesh_->geometric_model();
  Teuchos::RCP<const AmanziGeometry::Region> reg_ptr = gm_ptr->FindRegion(region_);

  if (reg_ptr->get_type() == AmanziGeometry::RegionType::POLYGON) {
    Teuchos::RCP<const AmanziGeometry::RegionPolygon> poly_reg =
      Teuchos::rcp_static_cast<const AmanziGeometry::RegionPolygon>(reg_ptr);
    reg_normal_ = poly_reg->normal();
    obs_planar_ = true;
  } else if (reg_ptr->get_type() == AmanziGeometry::RegionType::PLANE) {
    Teuchos::RCP<const AmanziGeometry::RegionPlane> plane_reg =
      Teuchos::rcp_static_cast<const AmanziGeometry::RegionPlane>(reg_ptr);
    reg_normal_ = plane_reg->normal();
    obs_planar_ = true;
  }

  if (variable_ == "aqueous mass flow rate" || variable_ == "aqueous volumetric flow rate" ||
      variable_ == "fractures aqueous volumetric flow rate") { // flux needs faces
    region_size_ = mesh_->get_set_size(
      region_, Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::OWNED);
    entity_ids_.resize(region_size_);
    mesh_->get_set_entities_and_vofs(
      region_, AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED, &entity_ids_, &vofs_);
    obs_boundary_ = 1;
    for (int i = 0; i != region_size_; ++i) {
      int f = entity_ids_[i];
      Amanzi::AmanziMesh::Entity_ID_List cells;
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      if (cells.size() == 2) {
        obs_boundary_ = 0;
        break;
      }
    }
    // to enforce common data on all processors
    int dummy(obs_boundary_);
    mesh_->get_comm()->MinAll(&dummy, &obs_boundary_, 1);

  } else { // all others need cells
    region_size_ = mesh_->get_set_size(region_, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    entity_ids_.resize(region_size_);
    mesh_->get_set_entities_and_vofs(
      region_, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &entity_ids_, &vofs_);
  }

  // find global mesh block size
  int dummy(region_size_);
  int global_mesh_block_size(0);
  mesh_->get_comm()->SumAll(&dummy, &global_mesh_block_size, 1);

  return global_mesh_block_size;
}


/* ******************************************************************
* Computes aqueous observations. Units should be taken from fields
* but fields do not populate them yet (FIXME).
****************************************************************** */
void
ObservableAqueous::ComputeObservation(State& S,
                                      double* value,
                                      double* volume,
                                      std::string& unit,
                                      double dt)
{
  Errors::Message msg;

  int dim = mesh_->space_dimension();

  // separate cases for density
  double rho = S.Get<double>("const_fluid_density");
  Key mol_density_key = Keys::getKey(domain_, "molar_density_liquid");

  Teuchos::RCP<const Epetra_MultiVector> rho_c;
  if (S.HasRecord(mol_density_key))
    rho_c = S.Get<CompositeVector>(mol_density_key).ViewComponent("cell");

  Key head_key = Keys::getKey(domain_, "hydraulic_head");
  Key poro_key = Keys::getKey(domain_, "porosity");
  Key sat_key = Keys::getKey(domain_, "saturation_liquid");
  Key wc_key = Keys::getKey(domain_, "water_content");
  Key pressure_key = Keys::getKey(domain_, "pressure");
  Key perm_key = Keys::getKey(domain_, "permeability");

  S.GetEvaluator(wc_key).Update(S, "cycle driver");
  const auto& wc = *S.Get<CompositeVector>(wc_key).ViewComponent("cell");
  const auto& porosity = *S.Get<CompositeVector>(poro_key).ViewComponent("cell");

  unit = "";

  if (variable_ == "volumetric water content") {
    for (int i = 0; i < region_size_; i++) {
      int c = entity_ids_[i];
      double vol = mesh_->cell_volume(c);
      *volume += vol;
      *value += wc[0][c] * vol;
    }
  } else if (variable_ == "gravimetric water content") {
    Key pd_key = Keys::getKey(domain_, "particle_density");
    if (!S.HasRecord(pd_key)) {
      msg << "Observation \"" << variable_ << "\" requires field \"particle_density\".\n";
      Exceptions::amanzi_throw(msg);
    }

    const auto& pd = *S.Get<CompositeVector>(pd_key).ViewComponent("cell");

    for (int i = 0; i < region_size_; i++) {
      int c = entity_ids_[i];
      double tmp = (rho_c.get()) ? (*rho_c)[0][c] / CommonDefs::MOLAR_MASS_H2O : rho;

      double vol = mesh_->cell_volume(c);
      *volume += vol;
      *value += wc[0][c] * tmp / (pd[0][c] * (1.0 - porosity[0][c])) * vol;
    }
  } else if (variable_ == "aqueous pressure") {
    const auto& pressure = *S.Get<CompositeVector>(pressure_key).ViewComponent("cell");

    for (int i = 0; i < region_size_; i++) {
      int c = entity_ids_[i];
      double vol = mesh_->cell_volume(c);
      *volume += vol;
      *value += pressure[0][c] * vol;
    }
    unit = "Pa";
  } else if (variable_ == "water table") {
    *value = CalculateWaterTable_(S, entity_ids_);
    *volume = 1.0;
    unit = "m";
  } else if (variable_ == "aqueous saturation") {
    const auto& ws = *S.Get<CompositeVector>(sat_key).ViewComponent("cell");

    for (int i = 0; i < region_size_; i++) {
      int c = entity_ids_[i];
      double vol = mesh_->cell_volume(c);
      *volume += vol;
      *value += ws[0][c] * vol;
    }
  } else if (variable_ == "ponded depth") {
    Key depth_key = Keys::getKey(domain_, "ponded_depth");
    const auto& depth = *S.Get<CompositeVector>(depth_key).ViewComponent("cell");

    for (int i = 0; i < region_size_; i++) {
      int c = entity_ids_[i];
      double vol = mesh_->cell_volume(c);
      *volume += vol;
      *value += depth[0][c] * vol;
    }
  } else if (variable_ == "velocity x") {
    Key vel_key = Keys::getKey(domain_, "velocity");
    const auto& vel = *S.Get<CompositeVector>(vel_key).ViewComponent("cell");

    for (int i = 0; i < region_size_; i++) {
      int c = entity_ids_[i];
      double vol = mesh_->cell_volume(c);
      *volume += vol;
      *value += vel[0][c] * vol;
    }
  } else if (variable_ == "hydraulic head") {
    const auto& hydraulic_head = *S.Get<CompositeVector>(head_key).ViewComponent("cell");

    for (int i = 0; i < region_size_; ++i) {
      int c = entity_ids_[i];
      double vol = mesh_->cell_volume(c);
      *volume += vol;
      *value += hydraulic_head[0][c] * vol;
    }
    unit = "m";
  } else if (variable_ == "permeability-weighted hydraulic head") {
    const auto& hydraulic_head = *S.Get<CompositeVector>(head_key).ViewComponent("cell");
    const auto& perm = *S.Get<CompositeVector>(perm_key).ViewComponent("cell");

    for (int i = 0; i < region_size_; ++i) {
      int c = entity_ids_[i];
      double vol = mesh_->cell_volume(c);
      double kxy = (dim == 2) ? perm[1][c] : std::pow(perm[1][c] * perm[2][c], 0.5);
      *volume += vol * kxy;
      *value += hydraulic_head[0][c] * vol * kxy;
    }
    unit = "m";
  } else if (variable_ == "drawdown") {
    const auto& hydraulic_head = *S.Get<CompositeVector>(head_key).ViewComponent("cell");

    for (int i = 0; i < region_size_; ++i) {
      int c = entity_ids_[i];
      double vol = mesh_->cell_volume(c);
      *volume += vol;
      *value += hydraulic_head[0][c] * vol;
    }
    unit = "m";

    // zero drawdown at time = t0 will be written directly to the file.
    // if (od.size() > 0) {
    //   *value = od.begin()->(*value) * (*volume) - (*value);
    // }
  } else if (variable_ == "permeability-weighted drawdown") {
    const auto& hydraulic_head = *S.Get<CompositeVector>(head_key).ViewComponent("cell");
    const auto& perm = *S.Get<CompositeVector>(perm_key).ViewComponent("cell");

    for (int i = 0; i < region_size_; ++i) {
      int c = entity_ids_[i];
      double vol = mesh_->cell_volume(c);
      double kxy = (dim == 2) ? perm[1][c] : std::pow(perm[1][c] * perm[2][c], 0.5);
      *volume += vol * kxy;
      *value += hydraulic_head[0][c] * vol * kxy;
    }
    unit = "m";

    // zero drawdown at time = t0 wil be written directly to the file.
    // if (od.size() > 0) {
    //   *value = od.begin()->(*value) * (*volume) - (*value);
    // }
  } else if (variable_ == "aqueous mass flow rate" || variable_ == "aqueous volumetric flow rate") {
    Key vol_flowrate_key = Keys::getKey(domain_, "volumetric_flow_rate");
    Teuchos::RCP<const Epetra_MultiVector> aperture_rcp;
    const auto& flowrate = *S.Get<CompositeVector>(vol_flowrate_key).ViewComponent("face");
    if (domain_ == "fracture")
      aperture_rcp = S.Get<CompositeVector>("fracture-aperture").ViewComponent("cell");
    const auto& fmap = *S.Get<CompositeVector>(vol_flowrate_key).Map().Map("face", true);
    Amanzi::AmanziMesh::Entity_ID_List cells;

    if (obs_boundary_ == 1) { // observation is on a boundary set
      for (int i = 0; i != region_size_; ++i) {
        int f = entity_ids_[i];
        mesh_->face_get_cells(f, Amanzi::AmanziMesh::Parallel_type::ALL, &cells);

        int sign, c = cells[0];
        mesh_->face_normal(f, false, c, &sign);
        double area = mesh_->face_area(f);
        double scale = 1.0;
        if (domain_ == "fracture") scale = (*aperture_rcp)[0][c];

        double tmp(1.0);
        if (variable_ == "aqueous mass flow rate")
          tmp = (rho_c.get()) ? (*rho_c)[0][c] / CommonDefs::MOLAR_MASS_H2O : rho;

        int g = fmap.FirstPointInElement(f);
        *value += sign * flowrate[0][g] * tmp * scale;
        *volume += area * scale;
      }
    } else if (obs_planar_) { // observation is on an interior planar set
      for (int i = 0; i != region_size_; ++i) {
        int f = entity_ids_[i];
        const AmanziGeometry::Point& face_normal = mesh_->face_normal(f);
        double area = mesh_->face_area(f);
        double sign = (reg_normal_ * face_normal) / area;

        mesh_->face_get_cells(f, Amanzi::AmanziMesh::Parallel_type::ALL, &cells);
        int c = cells[0];

        double scale = 1.0;
        if (domain_ == "fracture") scale = (*aperture_rcp)[0][c];

        double tmp(1.0);
        if (variable_ == "aqueous mass flow rate")
          tmp = (rho_c.get()) ? (*rho_c)[0][c] / CommonDefs::MOLAR_MASS_H2O : rho;

        int g = fmap.FirstPointInElement(f);
        *value += sign * flowrate[0][g] * tmp * scale;
        *volume += area * scale;
      }
    } else {
      msg << "Observations of \"aqueous mass flow rate\" and \"aqueous volumetric flow rate\""
          << " are only possible for Polygon, Plane and Boundary side sets";
      Exceptions::amanzi_throw(msg);
    }
    unit = "kg/s";

  } else if (variable_ == "pH") {
    Key ph_key = Keys::getKey(domain_, "pH");
    const auto& pH = *S.Get<CompositeVector>(ph_key).ViewComponent("cell");

    for (int i = 0; i < region_size_; ++i) {
      int c = entity_ids_[i];
      double vol = mesh_->cell_volume(c);
      *volume += vol;
      *value += pH[0][c] * vol;
    }
  } else {
    msg << "Cannot make an observation for aqueous variable \"" << variable_ << "\"";
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
 * Auxiliary routine: calculate maximum water table in a region.
 ****************************************************************** */
double
ObservableAqueous::CalculateWaterTable_(State& S, AmanziMesh::Entity_ID_List& ids)
{
  auto pressure = S.Get<CompositeVector>("pressure").ViewComponent("cell", true);
  double patm = S.Get<double>("atmospheric_pressure");

  // initilize and apply the reconstruction operator
  Teuchos::ParameterList plist;
  Operators::ReconstructionCellLinear lifting(mesh_);

  lifting.Init(plist);
  lifting.Compute(ids, pressure, 0);

  // set up extreme values for water table
  int dim = mesh_->space_dimension();
  double zmin(1e+99), zmax(-1e+99), pref(-1e+99), value(-1e+99);

  // estimate water table
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  int found(0);
  for (int i = 0; i < ids.size(); i++) {
    int c = ids[i];
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    double pc = (*pressure)[0][c];
    pref = pc;

    mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);
    for (int n = 0; n < faces.size(); ++n) {
      const AmanziGeometry::Point& xf = mesh_->face_centroid(faces[n]);
      zmin = std::min(zmin, xf[dim - 1]);
      zmax = std::max(zmax, xf[dim - 1]);

      double pf = lifting.getValue(c, xf);
      double dp = pf - pc;

      if ((pf - patm) * (pc - patm) <= 0.0) {
        if (fabs(dp) > 1e-8) {
          double a = (patm - pc) / dp;
          value = xc[dim - 1] + (xf[dim - 1] - xc[dim - 1]) * a;
          found = 1;
          break;
        }
      }
    }
  }

  // parallel update
  double tmp_loc[3] = { value, pref, zmax };
  double tmp_glb[3];
  mesh_->get_comm()->MaxAll(tmp_loc, tmp_glb, 3);
  value = tmp_glb[0];
  pref = tmp_glb[1];
  zmax = tmp_glb[2];

  double zmin_tmp(zmin);
  mesh_->get_comm()->MinAll(&zmin_tmp, &zmin, 1);

  int found_tmp = found;
  mesh_->get_comm()->MaxAll(&found_tmp, &found, 1);

  // process fully saturated and dry cases
  if (found == 0) {
    if (pref < patm) value = zmin;
    if (pref > patm) value = zmax;
  }

  return value;
}

} // namespace Amanzi
