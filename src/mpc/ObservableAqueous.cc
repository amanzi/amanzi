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

#include <map>
#include "ObservableAqueous.hh"
#include "RegionPlane.hh"
#include "RegionPolygon.hh"
#include "ReconstructionCell.hh"
#include "Units.hh"

namespace Amanzi{

/* ******************************************************************
* Delegating constructor
****************************************************************** */
ObservableAqueous::ObservableAqueous(std::string variable,
                                     std::string region,
                                     std::string functional,
                                     Teuchos::ParameterList& plist,
                                     Teuchos::ParameterList& units_plist,
                                     Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
    Observable(variable, region, functional, plist, units_plist, mesh)
{};


/* ******************************************************************
* Defines  
****************************************************************** */
int ObservableAqueous::ComputeRegionSize()
{
  //int mesh_block_size;

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

  if (variable_ == "aqueous mass flow rate" ||
      variable_ == "aqueous volumetric flow rate" ||
      variable_ == "fractures aqueous volumetric flow rate" ||
      variable_ == "pressure head face")    
    {  // flux needs faces
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
         
  // find global mesh block size
  int dummy = region_size_; 
  int global_mesh_block_size(0);
  mesh_->get_comm()->SumAll(&dummy, &global_mesh_block_size, 1);
      
  return global_mesh_block_size;
}


/* ******************************************************************
* Computes aqueous observations. Units should be taken from fields 
* but fields do not populate them yet (FIXME).
****************************************************************** */
void ObservableAqueous::ComputeObservation(
   State& S, double* value, double* volume, std::string& unit)
{
  //double volume, value;
  Errors::Message msg;
  int dim = mesh_->space_dimension();
  double rho = *S.GetScalarData("fluid_density");
  const Epetra_MultiVector& porosity = *S.GetFieldData("porosity")->ViewComponent("cell");    
  const Epetra_MultiVector& ws = *S.GetFieldData("saturation_liquid")->ViewComponent("cell");
  const Epetra_MultiVector& pressure = *S.GetFieldData("pressure")->ViewComponent("cell");

  
  unit = "";

  if (variable_ == "volumetric water content") {
    for (int i = 0; i < region_size_; i++) {
      int c = entity_ids_[i];
      double vol = mesh_->cell_volume(c);
      *volume += vol;
      *value  += porosity[0][c] * ws[0][c] * vol;
    }
  } else if (variable_ == "gravimetric water content") {
    if (!S.HasField("particle_density")) {
      msg << "Observation \""  << variable_ << "\" requires field \"particle_density\".\n";
      Exceptions::amanzi_throw(msg);
    }
    const Epetra_MultiVector& pd = *S.GetFieldData("particle_density")->ViewComponent("cell");    
    
    for (int i = 0; i < region_size_; i++) {
      int c = entity_ids_[i];
      double vol = mesh_->cell_volume(c);
      *volume += vol;
      *value  += porosity[0][c] * ws[0][c] * rho / (pd[0][c] * (1.0 - porosity[0][c])) * vol;
    }    
  } else if (variable_ == "aqueous pressure") {
    for (int i = 0; i < region_size_; i++) {
      int c = entity_ids_[i];
      double vol = mesh_->cell_volume(c);
      *volume += vol;
      *value  += pressure[0][c] * vol;
    } 
    unit = "Pa";
  } else if (variable_ == "water table") {
    *value = CalculateWaterTable_(S, entity_ids_);
    *volume = 1.0;
    unit = "m";
  } else if (variable_ == "pressure head face") {
    const Epetra_MultiVector& pressure = *S.GetFieldData("pressure")->ViewComponent("face");
    const Epetra_MultiVector& pressure_c = *S.GetFieldData("pressure")->ViewComponent("cell");    
    const auto& fmap = *S.GetFieldData("pressure")->Map().Map("face", true);
    Amanzi::AmanziMesh::Entity_ID_List cells;
    
    //double grav = 9.80;
    for (int i = 0; i < region_size_; i++) {
      int f = entity_ids_[i];
      mesh_->face_get_cells(f, Amanzi::AmanziMesh::Parallel_type::ALL, &cells);
      
      int c = cells[0];
      double sq = mesh_->face_area(f);
      int g = fmap.FirstPointInElement(f);        
      *volume += sq;
      *value  += pressure[0][g] * sq;
      //*value += pressure_c[0][c] * sq;
      
      //std::cout << region_<< ": "<<pressure[0][g] << " "<<pressure[0][g] /  (grav * rho) <<" "<<sq<<" "<<*volume<<"\n";
    }     
    unit = "m";
    
  } else if (variable_ == "aqueous saturation") {
    for (int i = 0; i < region_size_; i++) {
      int c = entity_ids_[i];
      double vol = mesh_->cell_volume(c);
      *volume += vol;
      *value  += ws[0][c] * vol;
    }    
  } else if (variable_ == "hydraulic head") {
    const Epetra_MultiVector& hydraulic_head = *S.GetFieldData("hydraulic_head")->ViewComponent("cell");
 
    for (int i = 0; i < region_size_; ++i) {
      int c = entity_ids_[i];
      double vol = mesh_->cell_volume(c);
      *volume += vol;
      *value  += hydraulic_head[0][c] * vol;
    }
    unit = "m";
  } else if (variable_ == "permeability-weighted hydraulic head") {
    const Epetra_MultiVector& hydraulic_head = *S.GetFieldData("hydraulic_head")->ViewComponent("cell");
    const Epetra_MultiVector& perm = *S.GetFieldData("permeability")->ViewComponent("cell");

    for (int i = 0; i < region_size_; ++i) {
      int c = entity_ids_[i];
      double vol = mesh_->cell_volume(c);
      double kxy = (dim == 2) ? perm[1][c] : std::pow(perm[1][c] * perm[2][c], 0.5);
      *volume += vol * kxy;
      *value  += hydraulic_head[0][c] * vol * kxy;
    }
    unit = "m";
  } else if (variable_ == "drawdown") {
    const Epetra_MultiVector& hydraulic_head = *S.GetFieldData("hydraulic_head")->ViewComponent("cell");

    for (int i = 0; i < region_size_; ++i) {
      int c = entity_ids_[i];
      double vol = mesh_->cell_volume(c);
      *volume += vol;
      *value  += hydraulic_head[0][c] * vol;
    }
    unit = "m";

    // zero drawdown at time = t0 will be written directly to the file.
    // if (od.size() > 0) { 
    //   *value = od.begin()->(*value) * (*volume) - (*value);
    // }
  } else if (variable_ == "permeability-weighted drawdown") {
    const Epetra_MultiVector& hydraulic_head = *S.GetFieldData("hydraulic_head")->ViewComponent("cell");
    const Epetra_MultiVector& perm = *S.GetFieldData("permeability")->ViewComponent("cell");

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
  } else if (variable_ == "aqueous mass flow rate" || 
             variable_ == "aqueous volumetric flow rate") {
    double density(1.0);
    if (variable_ == "aqueous mass flow rate") density = rho;
    const Epetra_MultiVector& darcy_flux = *S.GetFieldData("darcy_flux")->ViewComponent("face");
    const auto& fmap = *S.GetFieldData("darcy_flux")->Map().Map("face", true);
    
    if (obs_boundary_) { // observation is on a boundary set
      Amanzi::AmanziMesh::Entity_ID_List cells;

      for (int i = 0; i != region_size_; ++i) {
        int f = entity_ids_[i];
        mesh_->face_get_cells(f, Amanzi::AmanziMesh::Parallel_type::ALL, &cells);

        int sign, c = cells[0];
        const auto& normal = mesh_->face_normal(f, false, c, &sign);
        double area = mesh_->face_area(f);
        int g = fmap.FirstPointInElement(f);        
        *value  += sign * darcy_flux[0][g] * density;
        *volume += area;
        std::cout<< region_ <<" face "<<f<<": darcy "<<darcy_flux[0][g]<<" "<<sign<<" "<<density<<" "<<*volume<<"\n";
      }
    } else if (obs_planar_) {  // observation is on an interior planar set
      for (int i = 0; i != region_size_; ++i) {
        int f = entity_ids_[i];
        const AmanziGeometry::Point& face_normal = mesh_->face_normal(f);
        double area = mesh_->face_area(f);
        double sign = reg_normal_ * face_normal / area;
    
        *value  += sign * darcy_flux[0][f] * density;
        *volume += area;
      }
    } else {
      msg << "Observations of \"aqueous mass flow rate\" and \"aqueous volumetric flow rate\""
          << " are only possible for Polygon, Plane and Boundary side sets";
      Exceptions::amanzi_throw(msg);
    }
    unit = "kg/s";

  // fractures
  } else if (variable_ == "fractures aqueous volumetric flow rate") {
    const Epetra_MultiVector& darcy_flux = *S.GetFieldData("darcy_flux")->ViewComponent("face");
    const Epetra_MultiVector& aperture = *S.GetFieldData("fracture-aperture")->ViewComponent("cell");

    if (obs_boundary_) {
      Amanzi::AmanziMesh::Entity_ID_List cells;

      for (int i = 0; i != region_size_; ++i) {
        int f = entity_ids_[i];
        mesh_->face_get_cells(f, Amanzi::AmanziMesh::Parallel_type::ALL, &cells);

        int sign, c = cells[0];
        const auto& normal = mesh_->face_normal(f, false, c, &sign);
        double area = mesh_->face_area(f);

        *value  += sign * darcy_flux[0][f] * aperture[0][c];
        *volume += area * aperture[0][c];
      }
    } else {
      msg << "Observation \"" << variable_ << "\" is only possible for boundary side sets";
      Exceptions::amanzi_throw(msg);
    }
    unit = "kg/s";

  } else if (variable_ == "pH") {
    const Epetra_MultiVector& pH = *S.GetFieldData("pH")->ViewComponent("cell");

    for (int i = 0; i < region_size_; ++i) {
      int c = entity_ids_[i];
      double vol = mesh_->cell_volume(c);
      *volume += vol;
      *value  += pH[0][c] * vol;
    }
  } else {
    msg << "Cannot make an observation for aqueous variable \"" << variable_ << "\"";
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
 * Auxiliary routine: calculate maximum water table in a region.
 ****************************************************************** */
double ObservableAqueous::CalculateWaterTable_(State& S, 
                                               AmanziMesh::Entity_ID_List& ids)
{
  Teuchos::RCP<const Epetra_MultiVector> pressure = S.GetFieldData("pressure")->ViewComponent("cell", true);
  double patm = *S.GetScalarData("atmospheric_pressure");

  // initilize and apply the reconstruction operator
  Teuchos::ParameterList plist;
  Operators::ReconstructionCell lifting(mesh_);

  lifting.Init(pressure, plist);
  lifting.ComputeGradient(ids);
  const auto& gradient = *lifting.gradient(); 

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
  double tmp_loc[3] = {value, pref, zmax};
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

}  // namespace Amanzi

