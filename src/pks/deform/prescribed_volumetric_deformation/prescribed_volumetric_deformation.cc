/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Markus Berndt

   Interface for the Prescribed Deformation PK.

   <ParameterList name="prescribed deformation">
   <Parameter name="PK model" type="string" value="Prescibed Mesh Deformation"/>
   <Parameter name="Deformation method" type="string" value="method name"/>
   </ParameterList>

   ------------------------------------------------------------------------- */

#include "prescribed_volumetric_deformation.hh"
#include "porosity_evaluator.hh"

namespace Amanzi {
namespace Deform {

using namespace Amanzi::AmanziMesh;

RegisteredPKFactory<PrescribedVolumetricDeformation> PrescribedVolumetricDeformation::reg_("prescribed volumetric deformation");

PrescribedVolumetricDeformation::PrescribedVolumetricDeformation(const Teuchos::RCP<Teuchos::ParameterList>& plist,
        Teuchos::ParameterList& FElist,
        const Teuchos::RCP<TreeVector>& solution):
    PKDefaultBase(plist, FElist, solution),
    PKPhysicalBase(plist, FElist, solution)
{
  poro_key_ = plist.get<std::string>("porosity key","porosity");

  bottom_surface_ = plist_->get<Teuchos::Array<std::string> >("bottom surface").toVector();
  deform_region_ = plist_->get<std::string>("region");
  deform_factor_ = plist_->get<double>("factor");
  deform_time_start_ = plist_->get<double>("start time");
  deform_time_end_ = plist_->get<double>("end time");
}

// -- Setup data
void PrescribedVolumetricDeformation::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBase::setup(S);

  std::vector<AmanziMesh::Entity_kind> location(1);
  location[0] = AmanziMesh::CELL;
  std::vector<int> num_dofs(1);
  num_dofs[0] = 1;
  std::vector<std::string> name(1);
  name[0] = "cell";

  S->RequireField(key_, name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponents(name, location, num_dofs);

  // create storage for the centroid displacement
  num_dofs[0] = 3;

  S->RequireField("centroid", name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponents(name, location, num_dofs);


  // create storage for the vertex coordinates
  // we need to checkpoint those to be able to create
  // the deformed mesh after restart
  location[0] = AmanziMesh::NODE;
  num_dofs[0] = 3;
  name[0] = "node";
  S->RequireField("vertex coordinate", name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponents(name, location, num_dofs);

  // Also required are porosity and cell_volume fields.
  S->RequireFieldEvaluator("cell_volume");

  S->RequireFieldEvaluator("porosity");
  S->RequireField("porosity")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell",AmanziMesh::CELL,1);
}

// -- Initialize owned (dependent) variables.
void PrescribedVolumetricDeformation::initialize(const Teuchos::Ptr<State>& S) {
  PKPhysicalBase::initialize(S);

  // the PK's initial condition sets the initial porosity.  From this, we
  // calculate the actual initial condition, which is the rock volume.
  Epetra_MultiVector& deform = *S->GetFieldData(key_,name_)
      ->ViewComponent("cell",false);
  int ncells = deform.MyLength();
  for (int c=0; c!=ncells; ++c) {
    deform[0][c] = (1.0 - deform[0][c]) * mesh_->cell_volume(c);
  }

  // initialize the centroid displacement to be zero
  Epetra_MultiVector& centroid = *S->GetFieldData("centroid",name_)
      ->ViewComponent("cell",false);

  for (int c=0; c != ncells; ++c) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    centroid[0][c] = xc[0];
    centroid[1][c] = xc[1];
    centroid[2][c] = xc[2];
  }
  S->GetField("centroid",name_)->set_initialized();

  // initialize the vertex coordinate to the current mesh
  AmanziGeometry::Point_List pos;
  Entity_ID_List nodeids;
  Amanzi::AmanziGeometry::Point coords;

  // spatial dimension
  int dim = mesh_->space_dimension();
  // number of vertices
  int nV = mesh_->num_entities(Amanzi::AmanziMesh::NODE,
          Amanzi::AmanziMesh::OWNED);
  coords.init(dim);

  Epetra_MultiVector& vc = *S->GetFieldData("vertex coordinate",name_)
      ->ViewComponent("node",false);

  // search the id of the mid point on the top
  for (int iV=0; iV<nV; iV++) {
    // get the coords of the node
    mesh_->node_get_coordinates(iV,&coords);

    for ( int s=0; s<dim; ++s ) {
      vc[s][iV] = coords[s];
    }
  }
  S->GetField("vertex coordinate",name_)->set_initialized();
}

bool PrescribedVolumetricDeformation::advance(double dt) {
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Advancing deformation PK from time " << S_->time() << " to "
          << S_next_->time() << " with step size " << dt << std::endl;
    *out_ << "----------------------------------------------------------------" << std::endl;
  }


  AmanziMesh::Mesh * write_access_mesh_ =  const_cast<AmanziMesh::Mesh*>(&*S_next_->GetMesh());
  AmanziMesh::Mesh * write_access_surface_mesh_ =  const_cast<AmanziMesh::Mesh*>(&*S_next_->GetMesh("surface"));
  AmanziMesh::Mesh * write_access_surface3d_mesh_ =  const_cast<AmanziMesh::Mesh*>(&*S_next_->GetMesh("surface_3d"));


  double time1 = S_next_->time();
  double time0 = S_->time();


  AmanziGeometry::Point_List newpos, finpos;
  AmanziGeometry::Point_List surface_newpos, surface_finpos;
  AmanziGeometry::Point_List surface3d_newpos, surface3d_finpos;
  Entity_ID_List nodeids, surface_nodeids, surface3d_nodeids;
  Amanzi::AmanziGeometry::Point coords, new_coords;

  int nc = write_access_mesh_->num_entities(Amanzi::AmanziMesh::CELL,
                                            Amanzi::AmanziMesh::OWNED);
  
  std::vector<double> target_cell_volumes(nc,0.0);
  std::vector<double> min_cell_volumes(nc);

  
  // get the cell volumes from the previous time step
  const Epetra_MultiVector& cv = *S_->GetFieldData("cell_volume")
      ->ViewComponent("cell",false);

  const Epetra_MultiVector& poro = *S_->GetFieldData(poro_key_)
    ->ViewComponent("cell",false);

  for (int ic = 0; ic != nc; ++ic) {
    min_cell_volumes[ic] = cv[0][ic]*(1.0 - poro[0][ic]);
  }


  const Epetra_MultiVector& deform = *S_->GetFieldData(name_)
    ->ViewComponent("cell",false);

  double fac(1.0);

  double T0(std::max(deform_time_start_, time0));
  double T1(std::min(deform_time_end_, time1));

  if (T0 < T1) {
    double fT0 = 1.0 - (1.0 - deform_factor_)/(deform_time_end_-deform_time_start_)*(T0-deform_time_start_);
    double fT1 = 1.0 - (1.0 - deform_factor_)/(deform_time_end_-deform_time_start_)*(T1-deform_time_start_);    
      
    fac = fT1 / fT0;
  }

  unsigned int mesh_block_size = write_access_mesh_->get_set_size(deform_region_,
                                                                  Amanzi::AmanziMesh::CELL,
                                                                  Amanzi::AmanziMesh::OWNED);
  Amanzi::AmanziMesh::Entity_ID_List cell_ids(mesh_block_size);
  write_access_mesh_->get_set_entities(deform_region_, Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED,
                                       &cell_ids);
  for( Amanzi::AmanziMesh::Entity_ID_List::iterator c = cell_ids.begin(); c != cell_ids.end();  c++) {
    target_cell_volumes[*c] = fac * cv[0][*c];
  }
  
  // deform the mesh
  write_access_mesh_->deform( target_cell_volumes, min_cell_volumes, bottom_surface_, true); // deforms the mesh itself

  // for( Amanzi::AmanziMesh::Entity_ID_List::iterator c = cell_ids.begin(); c != cell_ids.end();  c++) {
  //   std::cout << min_cell_volumes[*c] << " " << target_cell_volumes[*c] << " " << cv[0][*c] << " " << write_access_mesh_->cell_volume(*c) <<std::endl;
  // }


  // now we have to adapt the surface mesh to the new volume mesh
  // extract the correct new coordinates for the surface from the domain
  // mesh and update the surface mesh accordingly
  int nV = write_access_surface_mesh_->num_entities(Amanzi::AmanziMesh::NODE,
          Amanzi::AmanziMesh::OWNED);
  AmanziGeometry::Point coord_domain(3);
  AmanziGeometry::Point coord_surface(2);

  for (int iV=0; iV<nV; iV++) {
    // get the coords of the node
    AmanziGeometry::Point coord(3);
    int node = write_access_surface_mesh_->entity_get_parent(AmanziMesh::NODE, iV);
    write_access_mesh_->node_get_coordinates(node, &coord_domain);

    // surface points are two dimensional
    coord_surface[0] = coord_domain[0];
    coord_surface[1] = coord_domain[1];

    surface_nodeids.push_back(iV);
    surface_newpos.push_back(coord_surface);

    surface3d_nodeids.push_back(iV);
    surface3d_newpos.push_back(coord_domain);
  }
  // now deform the surface meshes, note that we set the keep_valid flag to false, since
  // we want the surface mesh to exactly mirror the domain mesh
  write_access_surface_mesh_->deform( surface_nodeids, surface_newpos, false, &surface_finpos);
  write_access_surface3d_mesh_->deform( surface3d_nodeids, surface3d_newpos, false, &surface_finpos);



  solution_evaluator_->SetFieldAsChanged(S_next_.ptr()); // mark the placeholder evaluator as changed 

  // update vertex coordinates in state
  Epetra_MultiVector& vc = *S_next_->GetFieldData("vertex coordinate",name_)
      ->ViewComponent("node",false);
  int dim = write_access_mesh_->space_dimension();
  nV = write_access_mesh_->num_entities(Amanzi::AmanziMesh::NODE,
          Amanzi::AmanziMesh::OWNED);
  // loop over vertices and update vc
  for (int iV=0; iV<nV; iV++) {
    // get the coords of the node
    write_access_mesh_->node_get_coordinates(iV,&coords);
    for ( int s=0; s<dim; ++s ) {
      vc[s][iV] = coords[s];
    }
  }

  // setting deformation to be rock_volume at old time, (1-poro_old)*CV_old
  Epetra_MultiVector& def = *S_next_->GetFieldData(key_,name_)
      ->ViewComponent("cell",false);

  int ncells = def.MyLength();
  for (int c=0; c!=ncells; ++c) {
    def[0][c] = (1. - poro[0][c]) * cv[0][c];
  }

  // store the current centroids
  Epetra_MultiVector& centroid = *S_next_->GetFieldData("centroid",name_)
      ->ViewComponent("cell",false);

  for (int c=0; c != ncells; ++c) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    centroid[0][c] = xc[0];
    centroid[1][c] = xc[1];
    centroid[2][c] = xc[2];
  }
  
  return false;
}




} // namespace
} // namespace
