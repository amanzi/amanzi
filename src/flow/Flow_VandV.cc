/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <set>

#include "errors.hh"
#include "Flow_PK.hh"


namespace Amanzi {
namespace AmanziFlow {

/* ****************************************************************
* Construct default state for unit tests. It completes
* initialization of all missed objects in the state.
**************************************************************** */
void Flow_PK::CreateDefaultState(Teuchos::RCP<const AmanziMesh::Mesh>& mesh) 
{
  if (!S_->HasField("fluid_density")) {
    S_->RequireScalar("fluid_density", passwd_);
  }
  if (!S_->HasField("fluid_voscosity")) {
    S_->RequireScalar("fluid_viscosity", passwd_);
  }
  if (!S_->HasField("gravity")) {
    S_->RequireConstantVector("gravity", passwd_, dim);
  }

  if (!S_->HasField("permeability")) {
    S_->RequireField("permeability", passwd_)->SetMesh(mesh)->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, dim);
  }
  if (!S_->HasField("porosity")) {
    S_->RequireField("porosity", passwd_)->SetMesh(mesh)->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
 
  if (!S_->HasField("water_saturation")) {
    S_->RequireField("water_saturation", passwd_)->SetMesh(mesh)->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S_->HasField("prev_water_saturation")) {
    S_->RequireField("prev_water_saturation", passwd_)->SetMesh(mesh_)->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S_->HasField("specific_storage")) {
    S_->RequireField("specific_storage", passwd_)->SetMesh(mesh_)->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  if (!S_->HasField("specific_yield")) {
    S_->RequireField("specific_yield", passwd_)->SetMesh(mesh_)->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  if (!S_->HasField("pressure")) {
    std::vector<std::string> names;  // pressure
    names.push_back("cell");
    names.push_back("face");

    std::vector<AmanziMesh::Entity_kind> locations(2);
    locations[0] = AmanziMesh::CELL;
    locations[1] = AmanziMesh::FACE;

    std::vector<int> dofs(2, 1);
    S_->RequireField("pressure", passwd_)->SetMesh(mesh_)->SetGhosted(true)
        ->SetComponents(names, locations, dofs);
  }
  if (!S_->HasField("hydraulic_head")) {
    S_->RequireField("hydraulic_head", passwd_)->SetMesh(mesh_)->SetGhosted(false)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }
  
  if (!S_->HasField("darcy_flux")) {
    S_->RequireField("darcy_flux", passwd_)->SetMesh(mesh_)->SetGhosted(true)
        ->SetComponent("face", AmanziMesh::FACE, 1);
  }
  if (!S_->HasField("darcy_velocity")) {
    S_->RequireField("darcy_velocity", passwd_)->SetMesh(mesh)->SetGhosted(true)
        ->SetComponent("cell", AmanziMesh::CELL, dim);
  }
  
  // initialize fields
  S_->Setup();

  // set popular default values
  S_->GetFieldData("porosity", passwd_)->PutScalar(0.2);
  S_->GetField("porosity", passwd_)->set_initialized();

  *(S_->GetScalarData("fluid_density", passwd_)) = 1000.0;
  S_->GetField("fluid_density", passwd_)->set_initialized();

  *(S_->GetScalarData("fluid_viscosity", passwd_)) = 0.001;
  S_->GetField("fluid_viscosity", passwd_)->set_initialized();

  Epetra_Vector& gvec = *S_->GetConstantVectorData("gravity", passwd_);
  gvec.PutScalar(0.0);
  gvec[dim - 1] = -9.80;
  S_->GetField("gravity", passwd_)->set_initialized();

  S_->GetFieldData("permeability", passwd_)->PutScalar(1.0);
  S_->GetField("permeability", passwd_)->set_initialized();

  S_->GetFieldData("water_saturation", passwd_)->PutScalar(1.0);
  S_->GetField("water_saturation", passwd_)->set_initialized();

  S_->GetFieldData("prev_water_saturation", passwd_)->PutScalar(1.0);
  S_->GetField("prev_water_saturation", passwd_)->set_initialized();

  S_->GetFieldData("specific_storage", passwd_)->PutScalar(0.0);
  S_->GetField("specific_storage", passwd_)->set_initialized();

  S_->GetFieldData("specific_yield", passwd_)->PutScalar(0.0);
  S_->GetField("specific_yield", passwd_)->set_initialized();

  S_->GetFieldData("pressure", passwd_)->PutScalar(0.0);
  S_->GetField("pressure", passwd_)->set_initialized();

  S_->GetFieldData("hydraulic_head", passwd_)->PutScalar(0.0);
  S_->GetField("hydraulic_head", passwd_)->set_initialized();

  S_->GetFieldData("darcy_flux", passwd_)->PutScalar(0.0);
  S_->GetField("darcy_flux", passwd_)->set_initialized();

  S_->GetFieldData("darcy_velocity", passwd_)->PutScalar(0.0);
  S_->GetField("darcy_velocity", passwd_)->set_initialized();

  S_->InitializeFields();
}


/* ******************************************************************
* TODO: Verify that a BC has been applied to every boundary face.
* Right now faces without BC are considered no-mass-flux.                                         
****************************************************************** */
void Flow_PK::ValidateBCs() const
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
 ****************************************************************** */
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


}  // namespace AmanziFlow
}  // namespace Amanzi

