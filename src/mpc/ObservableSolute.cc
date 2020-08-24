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
    
  if (reg_ptr->type() == AmanziGeometry::POLYGON) {
    Teuchos::RCP<const AmanziGeometry::RegionPolygon> poly_reg =
        Teuchos::rcp_static_cast<const AmanziGeometry::RegionPolygon>(reg_ptr);
    reg_normal_ = poly_reg->normal();
    obs_planar_ = true;
  } else if (reg_ptr->type() == AmanziGeometry::PLANE) {
    Teuchos::RCP<const AmanziGeometry::RegionPlane> plane_reg =
        Teuchos::rcp_static_cast<const AmanziGeometry::RegionPlane>(reg_ptr);
    reg_normal_ = plane_reg->normal();
    obs_planar_ = true;
  }

  std::string solute_var = comp_names_[tcc_index_] + " volumetric flow rate";

  if (variable_ == solute_var ||
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
    State& S, double* value, double* volume, std::string& unit)
{
  Errors::Message msg;

  Key ws_key =  Keys::getKey(domain_, "saturation_liquid");
  Key tcc_key = Keys::getKey(domain_, "total_component_concentration");
  Key poro_key = Keys::getKey(domain_, "porosity");
  Key darcy_key = Keys::getKey(domain_, "darcy_flux");
  
  // fields below are subject to various input conditions
  Key total_sorbed_key = Keys::getKey(domain_, "total_sorbed");

  if (!S.HasField(tcc_key)) {
    // bail out with default values if this field is not yet created
    *value = 0.0;
    *volume = 1.0;
    return;
  }

  const Epetra_MultiVector& ws = *S.GetFieldData(ws_key)->ViewComponent("cell");
  const Epetra_MultiVector& tcc = *S.GetFieldData(tcc_key)->ViewComponent("cell");
  const Epetra_MultiVector& porosity = *S.GetFieldData(poro_key)->ViewComponent("cell");    

  unit = units_.system().concentration;

  if (variable_ == comp_names_[tcc_index_] + " aqueous concentration") { 
    for (int i = 0; i < region_size_; i++) {
      int c = entity_ids_[i];
      double factor = porosity[0][c] * ws[0][c] * mesh_->cell_volume(c);
      factor *= units_.concentration_factor();

      *value += tcc[tcc_index_][c] * factor;
      *volume += factor;
    }

  } else if (variable_ == comp_names_[tcc_index_] + " sorbed concentration" &&
             S.HasField(total_sorbed_key)) {
    const auto& sorbed = *S.GetFieldData(total_sorbed_key)->ViewComponent("cell");

    // we assume constant particle density
    for (int i = 0; i < region_size_; i++) {
      int c = entity_ids_[i];
      double factor = (1.0 - porosity[0][c]) * mesh_->cell_volume(c);

      *value += sorbed[tcc_index_][c] * factor;
      *volume += factor;
    }

  } else if (variable_ == comp_names_[tcc_index_] + " free ion concentration") {
    Key free_ion_key = Keys::getKey(domain_, "primary_free_ion_concentration_" + comp_names_[tcc_index_]);
    if (!S.HasField(free_ion_key)) {
      msg << "Observation: state does not have field \"" << variable_ << "\"";
      Exceptions::amanzi_throw(msg);
    }
    const auto& free_ion = *S.GetFieldData(free_ion_key)->ViewComponent("cell");

    for (int i = 0; i < region_size_; i++) {
      int c = entity_ids_[i];
      double factor = mesh_->cell_volume(c);

      *value += free_ion[0][c] * factor;
      *volume += factor;
    }

  } else if (variable_ == comp_names_[tcc_index_] + " gaseous concentration") { 
    for (int i = 0; i < region_size_; i++) {
      int c = entity_ids_[i];
      double factor = porosity[0][c] * (1.0 - ws[0][c]) * mesh_->cell_volume(c);
      factor *= units_.concentration_factor();

      *value += tcc[tcc_index_][c] * factor;
      *volume += factor;
    }

  } else if (variable_ == comp_names_[tcc_index_] + " volumetric flow rate") {

    const Epetra_MultiVector& darcy_flux = *S.GetFieldData(darcy_key)->ViewComponent("face");
    const auto& fmap = *S.GetFieldData(darcy_key)->Map().Map("face", true);
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

        *value += std::max(0.0, sign * darcy_flux[0][g]) * tcc[tcc_index_][c] * factor;
        *volume += area * factor;
      }

    } else if (obs_planar_) {  // observation is on an interior planar set
      for (int i = 0; i != region_size_; ++i) {
        int f = entity_ids_[i];
        mesh_->face_get_cells(f, Amanzi::AmanziMesh::Parallel_type::ALL, &cells);

        int csign, c = cells[0];
        const AmanziGeometry::Point& face_normal = mesh_->face_normal(f, false, c, &csign);
        if (darcy_flux[0][f] * csign < 0) c = cells[1];

        double area = mesh_->face_area(f);
        double sign = (reg_normal_ * face_normal) * csign / area;
        double factor = units_.concentration_factor();
        int g = fmap.FirstPointInElement(f);
        
        *value += sign * darcy_flux[0][g] * tcc[tcc_index_][c] * factor;
        *volume += area * factor;
      }

    } else {
      msg << "Observations of \"SOLUTE volumetric flow rate\""
          << " is only possible for Polygon, Plane and Boundary side sets";
      Exceptions::amanzi_throw(msg);
    }
  } else {
    msg << "Cannot make an observation for solute variable \"" << variable_ << "\"";
    Exceptions::amanzi_throw(msg);
  }    
}
 
}  // namespace Amanzi

