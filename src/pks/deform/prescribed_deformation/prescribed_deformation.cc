
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
    PKPhysicalBase(plist,solution),
    prescribed_deformation_case_(1)
{
  prescribed_deformation_case_ = plist.get<int>("deformation function",1);
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
}

// -- Initialize owned (dependent) variables.
void PrescribedDeformation::initialize(const Teuchos::Ptr<State>& S) {
  PKPhysicalBase::initialize(S); 
  // initialize the centroid displacement to be zero
  
  Epetra_MultiVector& centroid = *S->GetFieldData("centroid",name_)
      ->ViewComponent("cell",false);

  int c_owned = centroid.MyLength();
  for (int c=0; c != c_owned; ++c) {
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

bool PrescribedDeformation::advance(double dt) {

  AmanziMesh::Mesh * write_access_mesh_ =  const_cast<AmanziMesh::Mesh*>(&*S_next_->GetMesh());
  
  double ss = S_next_->time();
  double ss0 = S_->time();

  std::cout << "Advancing from " << ss0 << " to " << ss << std::endl;

  AmanziGeometry::Point_List newpos, finpos;
  Entity_ID_List nodeids;  
  Amanzi::AmanziGeometry::Point coords, new_coords;


  switch (prescribed_deformation_case_) {
    
    case 1:  // this is sinusoidal over time deformation in the negative z-direction 
      {
        double thickness = 0.9 + 0.1*std::cos( ss / (365.25*24*60*60) * 2.0*M_PI );   
        double thickness0 = 0.9 + 0.1*std::cos( ss0 / (365.25*24*60*60) * 2.0*M_PI );
        std::cout << std::setprecision(16);
        std::cout << "  thicknesses: " << thickness << ", " << thickness0 << std::endl;
        
        double factor = thickness/thickness0;
        std::cout << "SETTING factor = " << factor << std::endl;
        
        // get space dimensions
        int dim = write_access_mesh_->space_dimension();
        
        // number of vertices
        int nV = write_access_mesh_->num_entities(Amanzi::AmanziMesh::NODE, 
                                                  Amanzi::AmanziMesh::OWNED);
        
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
      }
      break;
      
    case 2: // this is denting the the mesh with a gaussian that is centered at (x,y) = (0,0)
      {
        tmax_ = 3e7;
        sigma_ = .1;
        mag_ = 1.0;
        
        // get space dimensions
        int dim = write_access_mesh_->space_dimension();
        
        // number of vertices
        int nV = write_access_mesh_->num_entities(Amanzi::AmanziMesh::NODE, 
                                                  Amanzi::AmanziMesh::OWNED);
              
        coords.init(dim);
        new_coords.init(dim);
        
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
      
  }
  

  // get the cell volumes from the previous time step
  const Epetra_MultiVector& cv = *S_->GetFieldData("cell_volume")
    ->ViewComponent("cell",false);
  
  // deform the mesh
  write_access_mesh_->deform( nodeids, newpos, true, &finpos); // deforms the mesh itself
  solution_evaluator_->SetFieldAsChanged(S_next_.ptr()); // mark the placeholder evaluator as changed


  // update vertex coordinates in state
  Epetra_MultiVector& vc = *S_next_->GetFieldData("vertex coordinate",name_)
    ->ViewComponent("node",false);
  int dim = write_access_mesh_->space_dimension();
  int nV = write_access_mesh_->num_entities(Amanzi::AmanziMesh::NODE, 
					    Amanzi::AmanziMesh::OWNED);  
  // loop over vertices and update vc
  for (int iV=0; iV<nV; iV++) {
    // get the coords of the node
    write_access_mesh_->node_get_coordinates(iV,&coords);
    for ( int s=0; s<dim; ++s ) {
      vc[s][iV] = coords[s];
    }
  }

  // -- this is no longer necessary, but I'm leaving it here to note how it
  //    should have been if it were needed...
  // get the new cell volumes, after deformation
  // bool changed = S_next_->GetFieldEvaluator("cell_volume")
  //     ->HasFieldChanged(S_next_.ptr(), name_);  // update the cell volume if the
  //                                               // placeholder eval has changed
  // const CompositeVector& cv_new = *S_next_->GetFieldData("cell_volume");
  // -- END no longer needed code.
  
  // setting deformation to be rock_volume at old time, (1-poro_old)*CV_old
  Epetra_MultiVector& def = *S_next_->GetFieldData("deformation",name_)
    ->ViewComponent("cell",false);
  const Epetra_MultiVector& poro = *S_->GetFieldData("porosity")
    ->ViewComponent("cell",false);
  
  int c_owned = def.MyLength();
  for (int c=0; c!=c_owned; ++c) {
    def[0][c] = (1. - poro[0][c]) * cv[0][c];
  }
  
  
  // store the current centroids
  Epetra_MultiVector& centroid = *S_next_->GetFieldData("centroid",name_)
    ->ViewComponent("cell",false);
  
  for (int c=0; c != c_owned; ++c) {
    const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
    centroid[0][c] = xc[0];
    centroid[1][c] = xc[1];
    centroid[2][c] = xc[2];
  }
  
  return false;
}
  
double PrescribedDeformation::deformation_fn_2 (double x, double y, double t) {

  if (t < tmax_) {
    return t/tmax_ * mag_*exp(-sigma_*(x*x+y*y));
  } else {
    return mag_*exp(-sigma_*(x*x+y*y));
  }

}





} // namespace
} // namespace

