/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Process kernel for coupling of Flow_PK and ShallowWater_PK.
*/

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
    Amanzi::PK_MPCSubcycled(pk_tree, global_list, S, soln) {
}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double SurfaceSubsurface_PK::get_dt() {
  //double dt = Amanzi::PK_MPCSubcycled::get_dt();
  double dt = sub_pks_[master_]->get_dt();
  return dt;
}


// -----------------------------------------------------------------------------
// Initialize
// -----------------------------------------------------------------------------
void SurfaceSubsurface_PK::Initialize(const Teuchos::Ptr<State>& S)
{
  
  mesh_domain_ = S->GetMesh();
  mesh_surface_ = S->GetMesh("surface");
  
  int nsurf_nodes = mesh_surface_->num_entities(AmanziMesh::NODE, AmanziMesh::Parallel_type::OWNED);
  int nsurf_cells = mesh_surface_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  
  Epetra_MultiVector &B_n = *S->GetFieldData("surface-bathymetry", "state")->ViewComponent("node");
  Epetra_MultiVector &B_c = *S->GetFieldData("surface-bathymetry", "state")->ViewComponent("cell");
  
  // Bathymetry node values
  for (int n = 0; n < nsurf_nodes; ++n) {
    int  v = mesh_surface_->entity_get_parent(AmanziMesh::Entity_kind::NODE, n);
    
    AmanziGeometry::Point node_crd;
    mesh_domain_->node_get_coordinates(v, &node_crd); // coordinate of surface node in 3D mesh
    
//    std::cout<<"z value: "<<node_crd[2]<<std::endl;
    B_n[0][n] = node_crd[2]; // z value of the surface
  }
  
  // Bathymetry cell values (compute from node values later)
  for (int c = 0; c < nsurf_cells; ++c) {
    int  cc = mesh_surface_->entity_get_parent(AmanziMesh::Entity_kind::CELL, c);

    const Amanzi::AmanziGeometry::Point &xc = mesh_domain_->cell_centroid(cc); // coordinate of surface node centroid in 3D mesh
    
//    std::cout<<"z cell value: "<<xc[2]<<std::endl;
    B_c[0][c] = xc[2]; // z value of the surface
  }
  
  S->GetField("surface-bathymetry", "state")->set_initialized();
  
  PK_MPC<PK>::Initialize(S);
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
void SurfaceSubsurface_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {
  sub_pks_[slave_]->CommitStep(t_old, t_new, S);
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool SurfaceSubsurface_PK::AdvanceStep(double t_old, double t_new, bool reinit) {
  bool fail = false;

  // advance the master PK using the full step size
  fail = sub_pks_[master_]->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;

  master_dt_ = t_new - t_old;

  sub_pks_[master_]->CommitStep(t_old, t_new, S_);

  slave_dt_ = sub_pks_[slave_]->get_dt();

  if (slave_dt_ > master_dt_) slave_dt_ = master_dt_;

  // advance the slave, subcycling if needed
  S_->set_intermediate_time(t_old);
  bool done = false;

  double dt_next = slave_dt_;
  double dt_done = 0.;
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
      // -- etc: unclear if state should be commited or not?
      // set the intermediate time
      S_->set_intermediate_time(t_old + dt_done + dt_next);
      sub_pks_[slave_]->CommitStep(t_old + dt_done, t_old + dt_done + dt_next, S_);
      dt_done += dt_next;
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
  return false;
}

}  // namespace Amanzi

