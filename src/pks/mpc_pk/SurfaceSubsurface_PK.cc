/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Process kernel for coupling of Flow_PK and ShallowWater_PK.
*/

#include "ShallowWater_PK.hh"
#include "SurfaceSubsurface_PK.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
SurfaceSubsurface_PK::SurfaceSubsurface_PK(
    Teuchos::ParameterList& pk_tree,
    const Teuchos::RCP<Teuchos::ParameterList>& global_list,
    const Teuchos::RCP<State>& S,
    const Teuchos::RCP<TreeVector>& soln) :
    Amanzi::PK_MPCSubcycled(pk_tree, global_list, S, soln)
{
  mesh_domain_ = S->GetMesh();
  mesh_surface_ = S->GetMesh("surface");
  
  Teuchos::ParameterList vlist;
  auto pk_list = Teuchos::sublist(global_list, "PKs", true)->sublist(name_);
  vlist.sublist("verbose object") = pk_list.sublist("verbose object");
  vo_ = Teuchos::rcp(new VerboseObject("SurfaceSubsurface", vlist)); 
}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double SurfaceSubsurface_PK::get_dt()
{
  // double dt = Amanzi::PK_MPCSubcycled::get_dt();
  double dt = sub_pks_[master_]->get_dt();
  return dt;
}


// -----------------------------------------------------------------------------
// Initialize
// -----------------------------------------------------------------------------
void SurfaceSubsurface_PK::Initialize()
{
  int nsurf_nodes = mesh_surface_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
  int nsurf_cells = mesh_surface_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  
  std::string passwd("");
  auto& B_n = *S_->GetW<CompositeVector>("surface-bathymetry", passwd).ViewComponent("node");
  auto& B_c = *S_->GetW<CompositeVector>("surface-bathymetry", passwd).ViewComponent("cell");
  
  // Bathymetry node values
  AmanziGeometry::Point node_crd;

  for (int n = 0; n < nsurf_nodes; ++n) {
    int v = mesh_surface_->entity_get_parent(AmanziMesh::Entity_kind::NODE, n);
    
    mesh_domain_->node_get_coordinates(v, &node_crd); // coordinate of surface node in 3D mesh
    
    B_n[0][n] = node_crd[2]; // z value of the surface
  }
  
  // Bathymetry cell values (compute from face centroids of the parent mesh)
  for (int c = 0; c < nsurf_cells; ++c) {
    int f = mesh_surface_->entity_get_parent(AmanziMesh::Entity_kind::CELL, c);

    const AmanziGeometry::Point &xf = mesh_domain_->face_centroid(f);
    
    B_c[0][c] = xf[2]; // z value of the surface
  }
  
  S_->GetRecordW("surface-bathymetry", passwd).set_initialized();
  
  PK_MPC<PK>::Initialize();
}


// -----------------------------------------------------------------------------
// Set master dt
// -----------------------------------------------------------------------------
void SurfaceSubsurface_PK::set_dt(double dt) {
  master_dt_ = dt;
  sub_pks_[master_]->set_dt(dt);
}


// -----------------------------------------------------------------------------
// Make necessary operatios by the end of the time steps.
// -----------------------------------------------------------------------------
void SurfaceSubsurface_PK::CommitStep(double t_old, double t_new, const Tag& tag) {
  sub_pks_[slave_]->CommitStep(t_old, t_new, tag);
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool SurfaceSubsurface_PK::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = false;

  // advance the master PK using the full step size
  fail = sub_pks_[master_]->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;

  master_dt_ = t_new - t_old;

  sub_pks_[master_]->CommitStep(t_old, t_new, Tags::DEFAULT);

  slave_dt_ = sub_pks_[slave_]->get_dt();

  if (slave_dt_ > master_dt_) slave_dt_ = master_dt_;

  // advance the slave, subcycling if needed
  S_->set_intermediate_time(t_old);

  int nsubcycles(0);
  bool done = false;
  double dt_next = slave_dt_;
  double dt_done(0.0), water_exchange(0.0);

  while (!done) {
    // do not overstep
    if (t_old + dt_done + dt_next > t_new) {
      dt_next = t_new - t_old - dt_done;
    }

    // take the step
    fail = sub_pks_[slave_]->AdvanceStep(t_old + dt_done, t_old + dt_done + dt_next, reinit);

    if (fail) {
      // if fail, cut the step and try again
      dt_next /= 2;
    } else {
      // if success, commit the state and increment to next intermediate
      S_->set_intermediate_time(t_old + dt_done + dt_next);
      sub_pks_[slave_]->CommitStep(t_old + dt_done, t_old + dt_done + dt_next, Tags::DEFAULT);

      dt_done += dt_next;
      nsubcycles++;

      // collect statistics
      auto sw_pk = Teuchos::rcp_dynamic_cast<ShallowWater::ShallowWater_PK>(sub_pks_[slave_]);
      water_exchange += sw_pk->get_total_source();

      // allow dt to grow only when success
      dt_next = sub_pks_[slave_]->get_dt();
    }

    // dt_next = sub_pks_[slave_]->get_dt();
    // no state recovery (e.g. pressure) is made, so the only option is to fail.
    if (dt_next < min_dt_)
       Exceptions::amanzi_throw("Failure in SurfaceSubsurface_PK: small time step.");

    // check for subcycling condition
    done = std::abs(t_old + dt_done - t_new) / (t_new - t_old) < 0.1*min_dt_;
  }

  // we reach this point when subcycling has been completed
  if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "shallow water PK made " << nsubcycles << " subcycles,  water exchange=" 
               << water_exchange << " m^3" << std::endl;
  }

  return false;
}

}  // namespace Amanzi

