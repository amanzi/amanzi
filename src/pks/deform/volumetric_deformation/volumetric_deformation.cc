/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Markus Berndt

   Interface for the Volumetric Deformation PK.

   <ParameterList name="volumetric deformation">
   <Parameter name="PK model" type="string" value="Prescibed Mesh Deformation"/>
   <Parameter name="Deformation method" type="string" value="method name"/>
   </ParameterList>

   ------------------------------------------------------------------------- */
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "composite_vector_function_factory.hh"

#include "volumetric_deformation.hh"
#include "porosity_evaluator.hh"

namespace Amanzi {
namespace Deform {

using namespace Amanzi::AmanziMesh;

RegisteredPKFactory<VolumetricDeformation> VolumetricDeformation::reg_("volumetric deformation");

VolumetricDeformation::VolumetricDeformation(Teuchos::ParameterList& plist,
        const Teuchos::RCP<TreeVector>& solution):
    PKDefaultBase(plist,solution),
    PKPhysicalBase(plist,solution) {
  poro_key_ = plist_.get<std::string>("porosity key","porosity");
  bottom_surface_name_ = plist_.get<std::string>("bottom surface region", "bottom face");
}

// -- Setup data
void VolumetricDeformation::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBase::setup(S);
  mesh_->build_columns(bottom_surface_name_);

  // save the meshes
  surf_mesh_ = S->GetMesh("surface");
  surf3d_mesh_ = S->GetMesh("surface_3d");
  mesh_nc_ = Teuchos::rcp_const_cast<AmanziMesh::Mesh>(mesh_);
  surf_mesh_nc_ = Teuchos::rcp_const_cast<AmanziMesh::Mesh>(surf_mesh_);
  surf3d_mesh_nc_ = Teuchos::rcp_const_cast<AmanziMesh::Mesh>(surf3d_mesh_);

  // create storage for primary variable, rock volume
  Teuchos::RCP<CompositeVectorFactory> cell_def_fac = S->RequireField(key_, name_);
  cell_def_fac->SetMesh(mesh_)->SetGhosted()
    ->SetComponent("cell", AmanziMesh::CELL, 1);

  // create the function
  Teuchos::ParameterList func_plist = plist_.sublist("deformation function");
  deform_func_ = Functions::CreateCompositeVectorFunction(func_plist, *cell_def_fac);

  // Additional data:
  // create storage for cell volume change
  S->RequireField("dcell_volume_dtime", name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("cell", AmanziMesh::CELL, 1);

  // create storage for the nodal deformation
  S->RequireField("nodal_dz", name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("node", AmanziMesh::NODE, 2);

  // create storage for the vertex coordinates
  // we need to checkpoint those to be able to create
  // the deformed mesh after restart
  int dim = mesh_->space_dimension();
  S->RequireField("vertex coordinate", name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("node", AmanziMesh::NODE, dim);

  S->RequireFieldEvaluator("cell_volume");
  S->RequireFieldEvaluator("porosity");
  S->RequireField("porosity")->SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell",AmanziMesh::CELL,1);
}

// -- Initialize owned (dependent) variables.
void VolumetricDeformation::initialize(const Teuchos::Ptr<State>& S) {
  PKPhysicalBase::initialize(S);

  // the PK's initial condition sets the initial porosity.  From this, we
  // calculate the actual initial condition, which is the rock volume.
  Epetra_MultiVector& rock_vol = *S->GetFieldData(key_,name_)
      ->ViewComponent("cell",false);
  int ncells = rock_vol.MyLength();
  for (int c=0; c!=ncells; ++c) {
    rock_vol[0][c] = (1.0 - rock_vol[0][c]) * mesh_->cell_volume(c);
  }

  // initialize the deformation rate to 0
  S->GetFieldData("dcell_volume_dtime",name_)->PutScalar(0.);
  S->GetField("dcell_volume_dtime",name_)->set_initialized();

  // initialize the plane displacement to be zero
  S->GetFieldData("nodal_dz",name_)->PutScalar(0.);
  S->GetField("nodal_dz",name_)->set_initialized();

  // initialize the vertex coordinate to the current mesh
  int dim = mesh_->space_dimension();
  AmanziGeometry::Point coords;
  coords.init(dim);

  // number of vertices
  int nV = mesh_->num_entities(Amanzi::AmanziMesh::NODE,
          Amanzi::AmanziMesh::OWNED);

  Epetra_MultiVector& vc = *S->GetFieldData("vertex coordinate",name_)
      ->ViewComponent("node",false);

  // search the id of the mid point on the top
  for (int iV=0; iV!=nV; ++iV) {
    // get the coords of the node
    mesh_->node_get_coordinates(iV,&coords);
    for (int s=0; s!=dim; ++s) vc[s][iV] = coords[s];
  }
  S->GetField("vertex coordinate",name_)->set_initialized();

}

bool VolumetricDeformation::advance(double dt) {
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_MEDIUM, true)) {
    *out_ << "Advancing deformation PK from time " << S_->time() << " to "
          << S_next_->time() << " with step size " << dt << std::endl;
    *out_ << "----------------------------------------------------------------" << std::endl;
  }

  double ss = S_next_->time();
  double ss0 = S_->time();
  double dT = ss - ss0;
  double T_mid = (ss + ss0) / 2.;

  // Evaluate the function at the mean time to get d(cv)/dT
  Teuchos::RCP<CompositeVector> dcell_vol_vec =
      S_next_->GetFieldData("dcell_volume_dtime", name_);
  deform_func_->Compute(T_mid, dcell_vol_vec.ptr());
  const Epetra_MultiVector& dcell_vol = *dcell_vol_vec->ViewComponent("cell",false);

  // Evaluate dz of the plane above the prism by integrating up the
  // change, and contribute these to the nodes that make up that plane
  // for averaging.
  Teuchos::RCP<CompositeVector> nodal_dz_vec = S_next_->GetFieldData("nodal_dz",name_);
  Epetra_MultiVector& nodal_dz_g = *nodal_dz_vec->ViewComponent("node",true);
  nodal_dz_g.PutScalar(0.);

  // Simultaneously evaluate contributions to the 

  AmanziMesh::Entity_ID_List bottomfaces;
  mesh_->get_set_entities(bottom_surface_name_, AmanziMesh::FACE,
                          AmanziMesh::OWNED, &bottomfaces);
  int ncols = 0;
  int ncells_tot = 0;

  for (unsigned int i=0; i!=bottomfaces.size(); ++i) {
    // find the bottom cell
    AmanziMesh::Entity_ID f_below = bottomfaces[i];
    AmanziMesh::Entity_ID_List cells;
    mesh_->face_get_cells(f_below, AmanziMesh::OWNED, &cells);
    ASSERT(cells.size() == 1);
    AmanziMesh::Entity_ID c0 = cells[0];

    // for each cell in the column, integrate the change
    double face_displacement = 0.;
    bool done = false;
    int ncells_in_col = 1;

    while (!done) {
      // find the face above
      AmanziMesh::Entity_ID f_above = mesh_->cell_get_face_above(c0);

      // calculate the displacement of the face above the cell
      double dz = mesh_->face_centroid(f_above)[2] - mesh_->face_centroid(f_below)[2];
      face_displacement += dz * dcell_vol[0][c0] * dT;

      // distribute that to the nodes of the face
      AmanziMesh::Entity_ID_List nodes;
      mesh_->face_get_nodes(f_above, &nodes);
      for (AmanziMesh::Entity_ID_List::const_iterator n=nodes.begin();
           n!=nodes.end(); ++n) {
        nodal_dz_g[0][*n] += face_displacement;
        nodal_dz_g[1][*n]++;
      }

      // iterate til we reach the top of the domain
      // get the cell above
      AmanziMesh::Entity_ID c1 = mesh_->cell_get_cell_above(c0);
      if (c1 >= 0) {
        c0 = c1;
        f_below = f_above;
        ncells_in_col++;
      } else {
        done = true;
      }
    }

    ncells_tot += ncells_in_col;
    if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
      *out_ << "  Integrated column " << ncols << " with " << ncells_in_col << " cells." << std::endl;
    }
    ncols++;
  }

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_EXTREME, true)) {
    int ncells_loc = ncells_tot;
    nodal_dz_vec->comm()->SumAll(&ncells_loc, &ncells_tot, 1);
    *out_ << "Integrated " << ncols << " columns, total # cells touched = " << ncells_tot << std::endl;
  }

  // take the average
  nodal_dz_vec->GatherGhostedToMaster("node");
  Epetra_MultiVector& nodal_dz_l = *nodal_dz_vec->ViewComponent("node",false);
  for (int i=0; i!=nodal_dz_l.MyLength(); ++i) {
    if (nodal_dz_l[1][i] > 1) {
      nodal_dz_l[0][i] /= nodal_dz_l[1][i];
    }
  }


  // form list of deformed nodes
  Entity_ID_List nodeids;
  AmanziGeometry::Point_List newpos;
  Amanzi::AmanziGeometry::Point coords, new_coords;
  int dim = mesh_->space_dimension();
  new_coords.init(dim);

  int nnodes = nodal_dz_l.MyLength();
  for (int i=0; i!=nnodes; ++i) {
    nodeids.push_back(i);
    mesh_->node_get_coordinates(i,&coords);
    new_coords = coords;
    new_coords[2] += nodal_dz_l[0][i];
    newpos.push_back(new_coords);
  }

  // deform the mesh
  AmanziGeometry::Point_List finpos;
  mesh_nc_->deform(nodeids, newpos, false, &finpos);

  // now we have to adapt the surface mesh to the new volume mesh
  // extract the correct new coordinates for the surface from the domain
  // mesh and update the surface mesh accordingly
  int nsurfnodes = surf_mesh_->num_entities(Amanzi::AmanziMesh::NODE,
          Amanzi::AmanziMesh::OWNED);
  AmanziGeometry::Point coord_domain(dim);
  AmanziGeometry::Point coord_surface(dim-1);

  Entity_ID_List surface_nodeids, surface3d_nodeids;
  AmanziGeometry::Point_List surface_newpos, surface3d_newpos;

  for (int i=0; i!=nsurfnodes; ++i) {
    // get the coords of the node
    AmanziGeometry::Point coord(3);
    AmanziMesh::Entity_ID pnode =
      surf_mesh_->entity_get_parent(AmanziMesh::NODE, i);
    mesh_->node_get_coordinates(i, &coord_domain);

    // surface points are two dimensional
    for (int s=0; s!=dim-1; ++s) coord_surface[s] = coord_domain[s];
    surface_nodeids.push_back(i);
    surface_newpos.push_back(coord_surface);

    surface3d_nodeids.push_back(i);
    surface3d_newpos.push_back(coord_domain);
  }

  // now deform the surface meshes, note that we set the keep_valid
  // flag to false, since we want the surface mesh to exactly mirror
  // the domain mesh
  AmanziGeometry::Point_List surface_finpos;
  surf_mesh_nc_->deform(surface_nodeids, surface_newpos, false, &surface_finpos);
  surf3d_mesh_nc_->deform(surface3d_nodeids, surface3d_newpos, false, &surface_finpos);

  // update vertex coordinates in state (for checkpointing)
  Epetra_MultiVector& vc = *S_next_->GetFieldData("vertex coordinate",name_)
      ->ViewComponent("node",false);
  for (int i=0; i!=nnodes; ++i) {
    mesh_->node_get_coordinates(i,&coords);
    for (int s=0; s!=dim; ++s) vc[s][i] = coords[s];
  }

  // setting deformation to be rock_volume at old time, (1-poro_old)*CV_old
  const Epetra_MultiVector& cv = *S_->GetFieldData("cell_volume")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& poro = *S_->GetFieldData(poro_key_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& def = *S_next_->GetFieldData(key_,name_)
      ->ViewComponent("cell",false);

  int ncells = def.MyLength();
  for (int c=0; c!=ncells; ++c) {
    def[0][c] = (1. - poro[0][c]) * cv[0][c];
  }
  // mark the placeholder evaluator as changed
  solution_evaluator_->SetFieldAsChanged(S_next_.ptr());

  // update porosity and cell volumes
  S_next_->GetFieldEvaluator("cell_volume")
      ->HasFieldChanged(S_next_.ptr(), name_);
  S_next_->GetFieldEvaluator("porosity")
      ->HasFieldChanged(S_next_.ptr(), name_);

  return false;
}


} // namespace
} // namespace
