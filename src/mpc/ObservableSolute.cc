/*
  Multi-Process Coordinator

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Konstantin Lipnikov
           Daniil Svyatsky
*/

#include "Evaluator.hh"

#include "ObservableSolute.hh"
#include "RegionPlane.hh"
#include "RegionPolygon.hh"

namespace Amanzi {

/* ******************************************************************
* Delegating constructor
****************************************************************** */
ObservableSolute::ObservableSolute(std::string variable,
                                   std::string region,
                                   std::string functional,
                                   Teuchos::ParameterList& plist,
                                   Teuchos::ParameterList& units_plist,
                                   Teuchos::RCP<const AmanziMesh::Mesh> mesh):
    Observable(variable, region, functional, plist, units_plist, mesh) 
{};


/* ******************************************************************
* TBW
****************************************************************** */
int ObservableSolute::ComputeRegionSize()
{
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

  if (variable_ == comp_names_[tcc_index_] + " volumetric flow rate" ||
      variable_ == comp_names_[tcc_index_] + " breakthrough curve" ||
      variable_ == "aqueous mass flow rate" ||
      variable_ == "aqueous volumetric flow rate") {  // flux needs faces
    region_size_ = mesh_->get_set_size(region_,
                                       Amanzi::AmanziMesh::FACE,
                                       Amanzi::AmanziMesh::Parallel_type::OWNED);
    entity_ids_.resize(region_size_);
    mesh_->get_set_entities_and_vofs(region_,
                                     AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED,
                                     &entity_ids_, 
                                     &vofs_);
    obs_boundary_ = true;
    for (int i = 0; i != region_size_; ++i) {
      int f = entity_ids_[i];
      Amanzi::AmanziMesh::Entity_ID_List cells;
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      if (cells.size() == 2) {
        obs_boundary_ = false;
        break;
      }
    }
  } else { // all others need cells
    region_size_ = mesh_->get_set_size(region_,
                                       AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);    
    entity_ids_.resize(region_size_);
    mesh_->get_set_entities_and_vofs(region_,
                                     AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED,
                                     &entity_ids_, &vofs_);
  }
      
  // find global meshblocksize
  int dummy = region_size_; 
  int global_mesh_block_size(0);
  mesh_->get_comm()->SumAll(&dummy, &global_mesh_block_size, 1);
      
  return global_mesh_block_size;
}


/* ******************************************************************
* Data calculation
****************************************************************** */
void ObservableSolute::ComputeObservation(
    State& S, double* value, double* volume, std::string& unit, double dt)
{
  Errors::Message msg;

  Key wc_key =  Keys::getKey(domain_, "water_content");
  Key ws_key =  Keys::getKey(domain_, "saturation_liquid");
  Key tcc_key = Keys::getKey(domain_, "total_component_concentration");
  Key poro_key = Keys::getKey(domain_, "porosity");
  Key darcy_key = Keys::getKey(domain_, "volumetric_flow_rate");
  
  // fields below are subject to various input conditions
  Key total_sorbed_key = Keys::getKey(domain_, "total_sorbed");

  if (!S.HasRecord(tcc_key)) {
    // bail out with default values if this field is not yet created
    *value = 0.0;
    *volume = 1.0;
    return;
  }

  const auto& tcc = *S.Get<CompositeVector>(tcc_key).ViewComponent("cell");
  const auto& porosity = *S.Get<CompositeVector>(poro_key).ViewComponent("cell");    

  unit = units_.system().concentration;

  if (variable_ == comp_names_[tcc_index_] + " aqueous concentration") { 
    S.GetEvaluator(wc_key).Update(S, "cycle driver");
    const auto& wc = *S.Get<CompositeVector>(wc_key).ViewComponent("cell");

    for (int i = 0; i < region_size_; i++) {
      int c = entity_ids_[i];
      double factor = wc[0][c] * mesh_->cell_volume(c);
      factor *= units_.concentration_factor();

      *value += tcc[tcc_index_][c] * factor;
      *volume += factor;
    }

  } else if (variable_ == comp_names_[tcc_index_] + " sorbed concentration" &&
             S.HasRecord(total_sorbed_key)) {
    const auto& sorbed = *S.Get<CompositeVector>(total_sorbed_key).ViewComponent("cell");

    // we assume constant particle density
    for (int i = 0; i < region_size_; i++) {
      int c = entity_ids_[i];
      double factor = (1.0 - porosity[0][c]) * mesh_->cell_volume(c);

      *value += sorbed[tcc_index_][c] * factor;
      *volume += factor;
    }

  } else if (variable_ == comp_names_[tcc_index_] + " free ion concentration") {
    Key free_ion_key = Keys::getKey(domain_, "primary_free_ion_concentration_" + comp_names_[tcc_index_]);
    if (!S.HasRecord(free_ion_key)) {
      msg << "Observation: state does not have field \"" << variable_ << "\"";
      Exceptions::amanzi_throw(msg);
    }
    const auto& free_ion = *S.Get<CompositeVector>(free_ion_key).ViewComponent("cell");

    for (int i = 0; i < region_size_; i++) {
      int c = entity_ids_[i];
      double factor = mesh_->cell_volume(c);

      *value += free_ion[0][c] * factor;
      *volume += factor;
    }

  } else if (variable_ == comp_names_[tcc_index_] + " gaseous concentration") {
    const auto& ws = *S.Get<CompositeVector>(ws_key).ViewComponent("cell");

    for (int i = 0; i < region_size_; i++) {
      int c = entity_ids_[i];
      double factor = porosity[0][c] * (1.0 - ws[0][c]) * mesh_->cell_volume(c);
      factor *= units_.concentration_factor();

      *value += tcc[tcc_index_][c] * factor;
      *volume += factor;
    }

  } else if (variable_ == comp_names_[tcc_index_] + " volumetric flow rate" ||
             variable_ == comp_names_[tcc_index_] + " breakthrough curve") {
    const auto& flowrate = *S.Get<CompositeVector>(darcy_key).ViewComponent("face");
    const auto& fmap = *S.Get<CompositeVector>(darcy_key).Map().Map("face", true);
    Amanzi::AmanziMesh::Entity_ID_List cells;

    if (obs_boundary_) { // observation is on a boundary set
      for (int i = 0; i != region_size_; ++i) {
        int f = entity_ids_[i];
        mesh_->face_get_cells(f, Amanzi::AmanziMesh::Parallel_type::ALL, &cells);

        int sign, c = cells[0];
        mesh_->face_normal(f, false, c, &sign);
        double area = mesh_->face_area(f);
        double factor = units_.concentration_factor();
        int g = fmap.FirstPointInElement(f);

        *value += std::max(0.0, sign * flowrate[0][g]) * tcc[tcc_index_][c] * factor;
        *volume += area * factor;
      }

    } else if (obs_planar_) {  // observation is on an interior planar set
      for (int i = 0; i != region_size_; ++i) {
        int f = entity_ids_[i];
        mesh_->face_get_cells(f, Amanzi::AmanziMesh::Parallel_type::ALL, &cells);

        int csign, c = cells[0];
        const AmanziGeometry::Point& face_normal = mesh_->face_normal(f, false, c, &csign);
        if (flowrate[0][f] * csign < 0) c = cells[1];

        double area = mesh_->face_area(f);
        double sign = (reg_normal_ * face_normal) * csign / area;
        double factor = units_.concentration_factor();
        int g = fmap.FirstPointInElement(f);
        
        *value += sign * flowrate[0][g] * tcc[tcc_index_][c] * factor;
        *volume += area * factor;
      }

    } else {
      msg << "Observations of \"SOLUTE volumetric flow rate\""
          << " is only possible for Polygon, Plane and Boundary side sets";
      Exceptions::amanzi_throw(msg);
    }

    if (variable_ == comp_names_[tcc_index_] + " breakthrough curve") {
      sum_ += *value * dt;
      *value = sum_;
    }
  } else {
    msg << "Cannot make an observation for solute variable \"" << variable_ << "\"";
    Exceptions::amanzi_throw(msg);
  }    
}
 
}  // namespace Amanzi

