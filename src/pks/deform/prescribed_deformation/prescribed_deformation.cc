/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
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

#include "prescribed_deformation.hh"
#include "porosity_evaluator.hh"

namespace Amanzi {
namespace Deform {

using namespace Amanzi::AmanziMesh;

PrescribedDeformation::PrescribedDeformation(const Teuchos::RCP<Teuchos::ParameterList>& plist,
        Teuchos::ParameterList& FElist,
        const Teuchos::RCP<TreeVector>& solution):
    PKDefaultBase(plist, FElist, solution),
    PKPhysicalBase(plist, FElist, solution),
    prescribed_deformation_case_(1)
{
  prescribed_deformation_case_ = plist_->get<int>("deformation function",1);
  poro_key_ = plist_->get<std::string>("porosity key","porosity");
}

// -- Setup data
void PrescribedDeformation::setup(const Teuchos::Ptr<State>& S) {
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
void PrescribedDeformation::initialize(const Teuchos::Ptr<State>& S) {
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

  // spatial dimension
  int dim = mesh_->space_dimension();
  Amanzi::AmanziGeometry::Point coords(dim);
  // number of vertices
  int nV = mesh_->num_entities(Amanzi::AmanziMesh::NODE,
          Amanzi::AmanziMesh::OWNED);

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

  // initialize functional parameters
  if (prescribed_deformation_case_ >= 2 &&
      prescribed_deformation_case_ <= 5) {
    tmax_ = plist_->get<double>("gaussian deformation time", 3.e7);
    sigma_ = plist_->get<double>("gaussian deformation sigma", 1.0);
    mag_ = plist_->get<double>("gaussian deformation amplitude", 3.0);
    z0_ = plist_->get<double>("initial height coordinate", 4.0);
    maxpert_ = plist_->get<double>("max magnitude perturbation", 0.0);
  }
}

bool PrescribedDeformation::advance(double dt) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "Advancing deformation PK from time " << S_->time() << " to "
               << S_next_->time() << " with step size " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;

  if (prescribed_deformation_case_ == 0) return false;

  AmanziMesh::Mesh * write_access_mesh_ =  const_cast<AmanziMesh::Mesh*>(&*S_next_->GetMesh());
  AmanziMesh::Mesh * write_access_surface_mesh_ =  const_cast<AmanziMesh::Mesh*>(&*S_next_->GetMesh("surface"));
  AmanziMesh::Mesh * write_access_surface3d_mesh_ =  const_cast<AmanziMesh::Mesh*>(&*S_next_->GetMesh("surface_3d"));


  double ss = S_next_->time();
  double ss0 = S_->time();

  //std::cout << "Advancing from " << ss0 << " to " << ss << std::endl;

  AmanziGeometry::Point_List newpos, finpos;
  AmanziGeometry::Point_List surface_newpos, surface_finpos;
  AmanziGeometry::Point_List surface3d_newpos, surface3d_finpos;
  Entity_ID_List nodeids, surface_nodeids, surface3d_nodeids;

  // get space dimensions
  int dim = write_access_mesh_->space_dimension();
  Amanzi::AmanziGeometry::Point coords(dim), new_coords(dim);

  switch (prescribed_deformation_case_) {

    case 1:  // this is sinusoidal over time deformation in the negative z-direction
      {
        double thickness = 0.9 + 0.1*std::cos( ss / (365.25*24*60*60) * 2.0*M_PI );
        double thickness0 = 0.9 + 0.1*std::cos( ss0 / (365.25*24*60*60) * 2.0*M_PI );
        if (vo_->os_OK(Teuchos::VERB_HIGH))
          *vo_->os() << std::setprecision(16)
                     << "  thicknesses: " << thickness << ", " << thickness0 << std::endl;

        double factor = thickness/thickness0;
        // std::cout << "SETTING factor = " << factor << std::endl;

        // number of vertices
        int nV = write_access_mesh_->num_entities(Amanzi::AmanziMesh::NODE,
                Amanzi::AmanziMesh::OWNED);

        // search the id of the mid point on the top
        for (int iV=0; iV<nV; iV++) {
          // get the coords of the node
          write_access_mesh_->node_get_coordinates(iV,&coords);

          for ( int s=0; s<dim-1; ++s ) {
            new_coords[s] = coords[s];
          }
          // new_coords[dim-1] = coords[dim-1] * fac;
          new_coords[dim-1] = coords[dim-1] * factor;

          // puch back for deform method
          nodeids.push_back(iV);
          newpos.push_back( new_coords );
        }
      }
      break;

    case 2: // this is denting the the mesh with a gaussian that is centered at (x,y) = (0,0)
      {
        // number of vertices
        int nV = write_access_mesh_->num_entities(Amanzi::AmanziMesh::NODE,
                Amanzi::AmanziMesh::OWNED);

        if (vo_->os_OK(Teuchos::VERB_HIGH))
          *vo_->os() << std::setprecision(16)
                     << "  z-coord (0,0) = " << deformation_fn_2(0.,0.,ss0) << std::endl;

        for (int iV=0; iV<nV; iV++) {
          // get the coords of the node
          write_access_mesh_->node_get_coordinates(iV,&coords);

          for ( int s=0; s<dim-1; ++s ) {
            new_coords[s] = coords[s];
          }

          // push the point down in z-direction
          new_coords[dim-1] = coords[dim-1] + deformation_fn_2(coords[0],coords[1],ss0) - deformation_fn_2(coords[0],coords[1],ss);

          // puch back for deform method
          nodeids.push_back(iV);
          newpos.push_back( new_coords );
        }
      }
      break;
  case 3: // this is denting the the mesh with a gaussian that is centered at (x,y) = (0,0)
  case 4:
  case 5:
      {
        // number of vertices
        int nV = write_access_mesh_->num_entities(Amanzi::AmanziMesh::NODE,
                Amanzi::AmanziMesh::OWNED);

        if (vo_->os_OK(Teuchos::VERB_HIGH))
          *vo_->os() << std::setprecision(16)
                     << "  z-coord (0,0) = " << deformation_fn_3(0.,0.,z0_,ss0) << std::endl;

        for (int iV=0; iV<nV; iV++) {
          // get the coords of the node
          write_access_mesh_->node_get_coordinates(iV,&coords);

          for ( int s=0; s<dim-1; ++s ) {
            new_coords[s] = coords[s];
          }

          // push the point down in z-direction
          if (prescribed_deformation_case_ == 3) {
            new_coords[dim-1] = coords[dim-1] + deformation_fn_3(coords[0],coords[1],coords[2],ss0) - deformation_fn_3(coords[0],coords[1],coords[2],ss);
          } else if (prescribed_deformation_case_ == 4) {
            new_coords[dim-1] = coords[dim-1] + deformation_fn_4(coords[0],coords[1],coords[2],ss0) - deformation_fn_4(coords[0],coords[1],coords[2],ss);
          } else if (prescribed_deformation_case_ == 5) {
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
  write_access_mesh_->deform( nodeids, newpos, true, &finpos); // deforms the mesh itself


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

double PrescribedDeformation::deformation_fn_2 (double x, double y, double t) {
  return std::min(1., t/tmax_) * mag_ * std::exp(- (x*x + y*y)/(2*sigma_*sigma_));
}

double PrescribedDeformation::deformation_fn_3 (double x, double y, double z, double t) {
  return std::min(1., t/tmax_) * mag_* z/z0_ * std::exp(- (x*x + y*y)/(2*sigma_*sigma_));
}

double PrescribedDeformation::deformation_fn_4 (double x, double y, double z, double t) {

  double fx0 = std::exp(- x*x / (2.0 *sigma_*sigma_) );
  double fx1 = std::exp(- (x-10.0)*(x-10.0) / (2.0 *sigma_*sigma_) );
  double fy0 = std::exp(- y*y / (2.0 *sigma_*sigma_) );
  double fy1 = std::exp(- (y-10.0)*(y-10.0) / (2.0 *sigma_*sigma_) );

  double fx0y0 = std::exp(- (x*x+y*y) / (2.0 *sigma_*sigma_) );
  double fx0y1 = std::exp(- (x*x+(y-10.0)*(y-10.0)) / (2.0 *sigma_*sigma_) ); 
  double fx1y1 = std::exp(- ((x-10.0)*(x-10.0)+(y-10.0)*(y-10.0)) / (2.0 *sigma_*sigma_) );  
  double fx1y0 = std::exp(- ((x-10.0)*(x-10.0)+y*y) / (2.0 *sigma_*sigma_) );    

  return std::min(1., t/tmax_) * mag_ * z/z0_ * (fx0+fx1+fy0+fy1-fx0y0-fx0y1-fx1y1-fx1y0);

}


double PrescribedDeformation::deformation_fn_5 (double x, double y, double z, double t) {

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
