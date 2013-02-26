
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

RegisteredPKFactory<PrescribedDeformation> PrescribedDeformation::reg_("prescribed deformation");

  
PrescribedDeformation::PrescribedDeformation(Teuchos::ParameterList& plist,
                                             const Teuchos::RCP<TreeVector>& solution):
    PKDefaultBase(plist,solution),
    PKPhysicalBase(plist,solution)
{



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

  location.resize(3);
  location[0] = AmanziMesh::CELL;  
  location[1] = AmanziMesh::CELL;
  location[2] = AmanziMesh::CELL;

  num_dofs.resize(3,1);
  
  name.resize(3);
  name[0] = "cell x";
  name[1] = "cell y";
  name[2] = "cell z";

  S->RequireField("centroid", name_)->SetMesh(mesh_)->SetGhosted()
    ->SetComponents(name, location, num_dofs);
  
}

// -- Initialize owned (dependent) variables.
void PrescribedDeformation::initialize(const Teuchos::Ptr<State>& S) {
  PKPhysicalBase::initialize(S); 
  // initialize the centroid displacement to be zero
  
  CompositeVector& centroid = *S->GetFieldData("centroid",name_);
  const AmanziMesh::Mesh& mesh = *S->GetMesh();

  int c_owned = centroid.size("cell x");
  for (int c=0; c != c_owned; ++c) {
    const AmanziGeometry::Point& xc = mesh.cell_centroid(c);
    centroid("cell x",c) = xc[0]; 
    centroid("cell y",c) = xc[1]; 
    centroid("cell z",c) = xc[2];
  }
  S->GetField("centroid",name_)->set_initialized();
}
  
bool PrescribedDeformation::advance(double dt) {

  AmanziMesh::Mesh * write_access_mesh_ =  const_cast<AmanziMesh::Mesh*>(&*S_next_->GetMesh());
  
  double ss = S_next_->time();
  double ss0 = S_->time();

  double thickness = 0.8 + 0.2*cos( ss / (365.25*24*60*60) * 2.0*M_PI );   
  double thickness0 = 0.8 + 0.2*cos( ss0 / (365.25*24*60*60) * 2.0*M_PI );

  double factor = thickness/thickness0;

  // get space dimensions
  int dim = write_access_mesh_->space_dimension();

  // set the list of the new position coords
  AmanziGeometry::Point_List newpos, finpos;
  
  // move the mid point on the mesh top
  Entity_ID_List nodeids;

  // number of vertices
  int nV = write_access_mesh_->num_entities(Amanzi::AmanziMesh::NODE, 
                               Amanzi::AmanziMesh::OWNED);

  Amanzi::AmanziGeometry::Point coords, new_coords;
  coords.init(dim);
  new_coords.init(dim);
  
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
  

  // get the cell volumes from the previous time step
  const CompositeVector& cv = *S_->GetFieldData("cell_volume"); 

  // deform the mesh
  write_access_mesh_->deform( nodeids, newpos, true, &finpos);

  // get the new cell volumes, after deformation
  const CompositeVector& cv_new = *S_next_->GetFieldData("cell_volume");

  // compute the cell-wise volume change factor
  CompositeVector& def = *S_next_->GetFieldData("deformation",name_);
  def.ReciprocalMultiply(1.0, cv, cv_new, 0.0);

  // store the current centroids
  CompositeVector& centroid = *S_next_->GetFieldData("centroid",name_);
  const AmanziMesh::Mesh& mesh = *S_next_->GetMesh();

  int c_owned = centroid.size("cell x");
  for (int c=0; c != c_owned; ++c) {
    const AmanziGeometry::Point& xc = mesh.cell_centroid(c);
    centroid("cell x",c) = xc[0]; 
    centroid("cell y",c) = xc[1]; 
    centroid("cell z",c) = xc[2];
  }  

  solution_evaluator_->SetFieldAsChanged();

}

} // namespace
} // namespace

