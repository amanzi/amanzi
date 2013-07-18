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
    PKPhysicalBase(plist,solution),
    volumetric_deformation_case_(1)
{
  poro_key_ = plist_.get<std::string>("porosity key","porosity");
  bottom_surface_name_ = plist_.get<std::string>("bottom surface region", "bottom face");
}

// -- Setup data
void VolumetricDeformation::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBase::setup(S);

  // create storage for primary variable, rock volume
  Teuchos::RCP<CompositeVectorFactory> cell_def_fac =
      S->RequireField(key_, name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("cell", AmanziMesh::CELL, 1);

  // create the function
  Teuchos::ParameterList func_plist = plist_.sublist("deformation function");
  deform_func_ = Functions::CreateCompositeVectorFunction(func_plist, cell_def_fac);

  // Additional data:
  // create storage for cell volume change
  S->RequireField("dcell_volume_dtime", name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("cell", AmanziMesh::CELL, 1);

  // create storage for the top plane's deformation
  S->RequireField("top_plane_dz", name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponents("cell", AmanziMesh::CELL, 1);

  // create storage for the vertex coordinates
  // we need to checkpoint those to be able to create
  // the deformed mesh after restart
  int dim = mesh_->space_dimension();
  S->RequireField("vertex coordinate", name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponents("node", AmanziMesh::NODE, dim);
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
  S->GetRecord("dcell_volume_dtime",name_)->set_initialized();

  // initialize the plane displacement to be zero
  S->GetFieldData("top_plane_dz",name_)->PutScalar(0.);
  S->GetRecord("top_plane_dz",name_)->set_initialized();

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

bool VolumetricDeformation::advance(double dt) {
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
    *out_ << "Advancing deformation PK from time " << S_->time() << " to "
          << S_next_->time() << " with step size " << dt << std::endl;
    *out_ << "----------------------------------------------------------------" << std::endl;
  }

  Teuchos::RCP<AmanziMesh::Mesh> write_access_mesh =
      Teuchos::rcp_const_cast<AmanziMesh::Mesh>(S_next_->GetMesh());
  Teuchos::RCP<AmanziMesh::Mesh> write_access_surface_mesh =
      Teuchos::rcp_const_cast<AmanziMesh::Mesh>(S_next_->GetMesh("surface"));
  Teuchos::RCP<AmanziMesh::Mesh> write_access_surface3d_mesh =
      Teuchos::rcp_const_cast<AmanziMesh::Mesh>(S_next_->GetMesh("surface_3d"));

  double ss = S_next_->time();
  double ss0 = S_->time();
  double dT = ss - ss0;
  double T_mid = (ss + ss0) / 2.;

  // Evaluate the function at the mean time to get d(cv)/dT
  Teuchos::RCP<CompositeVector> dcell_vol =
      S_next_->GetFieldData("dcell_volume_dtime", name_);
  deform_func_->Compute(T_mid, dcell_vol.ptr());

  // Evaluate dz of the plane above the prism by integrating up the change
  Epetra_MultiVector& plane_dz = *S_next_->GetFieldData("top_plane_dz",name_)
      ->ViewComponent("cell",false);

  AmanziMesh::Entity_ID_List bottomfaces;
  mesh_->get_set_entities(bottom_surface_name_, AmanziMesh::FACE,
                          AmanziMesh::OWNED, &bottomfaces);
  for (unsigned int i=0; i!=bottomfaces.size(); ++i) {
    // find the bottom cell
    AmanziMesh::Entity_ID f_below = bottomfaces[i];
    AmanziMesh::Entity_ID_List cells;
    mesh_->face_get_cells(f, AmanziMesh::OWNED, &cells);
    ASSERT(cells.size() == 1);
    AmanziMesh::Entity_ID c0 = cells[0];
    AmanziMesh::Entity_ID c1 = mesh_->cell_get_cell_above(c0);
    double deformation_below = 0.;

    // for each cell in the column, integrate the change
    while (c1 >= 0) {
      // find the common face
      std::vector<int> dirs;
      AmanziMesh::Entity_ID_List faces0 = mesh_->cell_get_faces_and_dirs(c0, &faces0, &dirs);
      AmanziMesh::Entity_ID_List faces1 = mesh_->cell_get_faces_and_dirs(c1, &faces1, &dirs);
      AmanziMesh::Entity_ID f_above = -1;

      for (AmanziMesh::Entity_ID_List::const_iterator f0=faces0.begin();
           f0!=faces0.end(); ++f0) {
        if (std::find(faces1.begin(), faces1.end(), *f0) - faces1.begin() < faces1.size()) {
          f_above = *f0;
          continue;
        }
      }
      ASSERT(common_face >= 0);

      // calculate the current height of the cell
      double cell_height = mesh_->face_centroid(f_above)[2] - mesh_->face_centroid(f_below)[2];

      // store the required deformation
      plane_dz[0][c0] = 

    }

  }


  AmanziGeometry::Point_List newpos, finpos;
  AmanziGeometry::Point_List surface_newpos, surface_finpos;
  AmanziGeometry::Point_List surface3d_newpos, surface3d_finpos;
  Entity_ID_List nodeids, surface_nodeids, surface3d_nodeids;
  Amanzi::AmanziGeometry::Point coords, new_coords;


  switch (volumetric_deformation_case_) {
    case 1: // calls function 1
      {
        // get space dimensions
        int dim = write_access_mesh->space_dimension();

        // number of vertices
        int nV = write_access_mesh->num_entities(Amanzi::AmanziMesh::NODE,
                Amanzi::AmanziMesh::OWNED);

        coords.init(dim);
        new_coords.init(dim);

        if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
          *out_ << std::setprecision(16);
          *out_ << "  z-coord (0,0) = " << deformation_fn_3(0.,0.,z0_,ss0) << std::endl;
        }

        for (int iV=0; iV<nV; iV++) {

          // get the coords of the node
          write_access_mesh->node_get_coordinates(iV,&coords);

          for ( int s=0; s<dim-1; ++s ) {
            new_coords[s] = coords[s];
          }

          // push the point down in z-direction
          if (volumetric_deformation_case_ == 3) {
            new_coords[dim-1] = coords[dim-1] + deformation_fn_3(coords[0],coords[1],coords[2],ss0) - deformation_fn_3(coords[0],coords[1],coords[2],ss);
          } else if (volumetric_deformation_case_ == 4) {
            new_coords[dim-1] = coords[dim-1] + deformation_fn_4(coords[0],coords[1],coords[2],ss0) - deformation_fn_4(coords[0],coords[1],coords[2],ss);
          } else if (volumetric_deformation_case_ == 5) {
            new_coords[dim-1] = coords[dim-1] + deformation_fn_5(coords[0],coords[1],coords[2],ss0) - deformation_fn_5(coords[0],coords[1],coords[2],ss);
          }

          // puch back for deform method
          nodeids.push_back(iV);
          newpos.push_back( new_coords );
        }
      }
      break;
  }


  // get the cell volumes from the previous time step
  const Epetra_MultiVector& cv = *S_->GetFieldData("cell_volume")
      ->ViewComponent("cell",false);

  // deform the mesh
  write_access_mesh->deform( nodeids, newpos, true, &finpos); // deforms the mesh itself


  // now we have to adapt the surface mesh to the new volume mesh
  // extract the correct new coordinates for the surface from the domain
  // mesh and update the surface mesh accordingly
  int nV = write_access_surface_mesh->num_entities(Amanzi::AmanziMesh::NODE,
          Amanzi::AmanziMesh::OWNED);
  AmanziGeometry::Point coord_domain(3);
  AmanziGeometry::Point coord_surface(2);

  for (int iV=0; iV<nV; iV++) {
    // get the coords of the node
    AmanziGeometry::Point coord(3);
    int node = write_access_surface_mesh->entity_get_parent(AmanziMesh::NODE, iV);
    write_access_mesh->node_get_coordinates(node, &coord_domain);

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
  write_access_surface_mesh->deform( surface_nodeids, surface_newpos, false, &surface_finpos);
  write_access_surface3d_mesh->deform( surface3d_nodeids, surface3d_newpos, false, &surface_finpos);



  solution_evaluator_->SetFieldAsChanged(S_next_.ptr()); // mark the placeholder evaluator as changed

  // update vertex coordinates in state
  Epetra_MultiVector& vc = *S_next_->GetFieldData("vertex coordinate",name_)
      ->ViewComponent("node",false);
  int dim = write_access_mesh->space_dimension();
  nV = write_access_mesh->num_entities(Amanzi::AmanziMesh::NODE,
          Amanzi::AmanziMesh::OWNED);
  // loop over vertices and update vc
  for (int iV=0; iV<nV; iV++) {
    // get the coords of the node
    write_access_mesh->node_get_coordinates(iV,&coords);
    for ( int s=0; s<dim; ++s ) {
      vc[s][iV] = coords[s];
    }
  }

  // setting deformation to be rock_volume at old time, (1-poro_old)*CV_old
  Epetra_MultiVector& def = *S_next_->GetFieldData(key_,name_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& poro = *S_->GetFieldData(poro_key_)
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


double VolumetricDeformation::deformation_fn_1(double x, double y, double z, double t) {

  double line1 = (0.5*x - y + 2.0)/sqrt(0.25+1);
  double line2 = (-2.0*x - y + 14.0)/sqrt(4.0+1.0);
  double fx0 = std::exp(- line1*line1 / (2.0 *sigma_*sigma_) );
  double fx1 = std::exp(- line2*line2 / (2.0 *sigma_*sigma_) );

  double xx = 24.0/5.0;
  double yy = 22.0/5.0;
  double fxy = std::exp(- ((x-xx)*(x-xx)+(y-yy)*(y-yy)) /(2.0 *sigma_*sigma_) );

  double varmag = mag_ *  ( 1.0 + maxpert_ * (x + y)/20.0 );

  return std::min(1., t/tmax_) * varmag * z/z0_ * (fx0+fx1 - fxy);
}



} // namespace
} // namespace
