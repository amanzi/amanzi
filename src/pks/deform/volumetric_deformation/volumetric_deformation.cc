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

#include "NKAOperator.hh"
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
}

// -- Setup data
void VolumetricDeformation::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBase::setup(S);

  // save the meshes
  mesh_nc_ = Teuchos::rcp_const_cast<AmanziMesh::Mesh>(mesh_);

  if (S->HasMesh("surface")) {
    surf_mesh_ = S->GetMesh("surface");
    surf3d_mesh_ = S->GetMesh("surface_3d");
    surf_mesh_nc_ = Teuchos::rcp_const_cast<AmanziMesh::Mesh>(surf_mesh_);
    surf3d_mesh_nc_ = Teuchos::rcp_const_cast<AmanziMesh::Mesh>(surf3d_mesh_);
  }

  // create storage for primary variable, rock volume
  Teuchos::RCP<CompositeVectorFactory> rv_fac =
      S->RequireField(key_, name_);
  rv_fac->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("cell", AmanziMesh::CELL, 1);

  // Create storage and a function for cell volume change
  Teuchos::RCP<CompositeVectorFactory> cv_fac =
      S->RequireField("dcell_volume_dtime", name_);
  cv_fac->Update(*rv_fac);
  Teuchos::ParameterList func_plist = plist_.sublist("deformation function");
  deform_func_ = Functions::CreateCompositeVectorFunction(func_plist, *cv_fac);

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

  // solve method
  Teuchos::ParameterList op_plist = plist_.sublist("global solve operator");

  // collect a set of the fixed nodes
  std::vector<std::string> bottom_regions;
  if (op_plist.isParameter("bottom region")) {
    bottom_regions.push_back(op_plist.get<std::string>("bottom region"));
  } else {
    bottom_regions =
        op_plist.get<Teuchos::Array<std::string> >("bottom regions").toVector();
  }
  std::string bottom_region_type = op_plist.get<std::string>("bottom region type","node");
  std::set<AmanziMesh::Entity_ID> bottom_node_set;
  for (std::vector<std::string>::const_iterator region=bottom_regions.begin();
       region!=bottom_regions.end(); ++region) {
    if (bottom_region_type == "node") {
      AmanziMesh::Entity_ID_List region_nodes;
      mesh_->get_set_entities(*region, AmanziMesh::NODE, AmanziMesh::OWNED,&region_nodes);
      bottom_node_set.insert(region_nodes.begin(), region_nodes.end());
    } else if (bottom_region_type == "face") {
      AmanziMesh::Entity_ID_List region_faces;
      mesh_->get_set_entities(*region, AmanziMesh::FACE, AmanziMesh::OWNED,&region_faces);
      for (AmanziMesh::Entity_ID_List::const_iterator f=region_faces.begin();
           f!=region_faces.end(); ++f) {
        AmanziMesh::Entity_ID_List face_nodes;
        mesh_->face_get_nodes(*f, &face_nodes);
        bottom_node_set.insert(face_nodes.begin(), face_nodes.end());
      }
    } else {
      Errors::Message mesg("Invalid bottom region type (must be node or face)");
      Exceptions::amanzi_throw(mesg);
    }
  }

  // create the unique list of owned bottom nodes
  unsigned int nnodes_owned = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::OWNED);
  Teuchos::RCP<AmanziMesh::Entity_ID_List> bottom_node_list =
      Teuchos::rcp(new AmanziMesh::Entity_ID_List());
  for (std::set<AmanziMesh::Entity_ID>::const_iterator n=bottom_node_set.begin();
       n!=bottom_node_set.end(); ++n) {
    if (*n < nnodes_owned) bottom_node_list->push_back(*n);
  }

  // create the operator
  def_matrix_ = Teuchos::rcp(new Operators::MatrixVolumetricDeformation(op_plist,
          mesh_, bottom_node_list));

  if (op_plist.isSublist("NKA Solver")) {
    Teuchos::ParameterList solver_plist = op_plist.sublist("NKA Solver");
    operator_ = Teuchos::rcp(
        new NKAOperator<CompositeMatrix,CompositeVector,
                        CompositeVectorFactory>(solver_plist,def_matrix_));
  } else {
    operator_ = def_matrix_;
  }

  // create storage for the nodal deformation
  S->RequireField("nodal_dz", name_)->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("node", AmanziMesh::NODE, 1);
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

  { // initialize the vertex coordinate to the current mesh
    int dim = mesh_->space_dimension();
    AmanziGeometry::Point coords;
    coords.init(dim);
    int nnodes = mesh_->num_entities(Amanzi::AmanziMesh::NODE,
            Amanzi::AmanziMesh::OWNED);

    Epetra_MultiVector& vc = *S->GetFieldData("vertex coordinate",name_)
        ->ViewComponent("node",false);

    // search the id of the mid point on the top
    for (int iV=0; iV!=nnodes; ++iV) {
      // get the coords of the node
      mesh_->node_get_coordinates(iV,&coords);
      for (int s=0; s!=dim; ++s) vc[s][iV] = coords[s];
    }
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

  const Epetra_MultiVector& cv = *S_->GetFieldData("cell_volume")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& poro = *S_->GetFieldData(poro_key_)
      ->ViewComponent("cell",false);

  // Evaluate the function at the mean time to get d(cv)/dT
  Teuchos::RCP<CompositeVector> dcell_vol_vec =
      S_next_->GetFieldData("dcell_volume_dtime", name_);
  dcell_vol_vec->PutScalar(0.);
  deform_func_->Compute(T_mid, dcell_vol_vec.ptr());

  { // scale by cell volume, this is fractional loss
    Epetra_MultiVector& dcell_vol_c = *dcell_vol_vec->ViewComponent("cell",false);
    for (int c=0; c!=dcell_vol_c.MyLength(); ++c) {
      dcell_vol_c[0][c] *= cv[0][c] * dT;
    }
  }

  // calculate the nodal deformation
  Teuchos::RCP<CompositeVector> nodal_dz_vec =
      S_next_->GetFieldData("nodal_dz",name_);
  Teuchos::RCP<CompositeVector> rhs = Teuchos::rcp(new CompositeVector(*nodal_dz_vec));
  def_matrix_->ApplyRHS(*dcell_vol_vec, rhs.ptr());
  operator_->ApplyInverse(*rhs, nodal_dz_vec.ptr());

  // form list of deformed nodes
  Entity_ID_List nodeids;
  AmanziGeometry::Point_List newpos;
  Amanzi::AmanziGeometry::Point coords, new_coords;
  int dim = mesh_->space_dimension();
  new_coords.init(dim);

  nodal_dz_vec->ScatterMasterToGhosted("node",true);

  // WORKAROUND for non-communication in deform() by Mesh
  //  const Epetra_MultiVector& nodal_dz_l =
  //      *nodal_dz_vec->ViewComponent("node",false);
  const Epetra_MultiVector& nodal_dz_l =
      *nodal_dz_vec->ViewComponent("node",true);

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

  if (surf_mesh_ != Teuchos::null) {
    // now we have to adapt the surface mesh to the new volume mesh
    // extract the correct new coordinates for the surface from the domain
    // mesh and update the surface mesh accordingly

    // WORKAROUND for non-communication in deform() by Mesh
    //    int nsurfnodes = surf_mesh_->num_entities(Amanzi::AmanziMesh::NODE,
    //            Amanzi::AmanziMesh::OWNED);
    int nsurfnodes = surf_mesh_->num_entities(Amanzi::AmanziMesh::NODE,
            Amanzi::AmanziMesh::USED);
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
  }

  {  // update vertex coordinates in state (for checkpointing)
    Epetra_MultiVector& vc = *S_next_->GetFieldData("vertex coordinate",name_)
        ->ViewComponent("node",false);
    int nnodes = vc.MyLength();
    for (int i=0; i!=nnodes; ++i) {
      mesh_->node_get_coordinates(i,&coords);
      for (int s=0; s!=dim; ++s) vc[s][i] = coords[s];
    }
  }

  // setting deformation to be rock_volume at old time, (1-poro_old)*CV_old
  {
    Epetra_MultiVector& def = *S_next_->GetFieldData(key_,name_)
        ->ViewComponent("cell",false);

    int ncells = def.MyLength();
    for (int c=0; c!=ncells; ++c) {
      def[0][c] = (1. - poro[0][c]) * cv[0][c];
    }
  }

  // mark the placeholder evaluator as changed
  solution_evaluator_->SetFieldAsChanged(S_next_.ptr());

  // REMOVE ME --etc
  // update porosity and cell volumes
  S_next_->GetFieldEvaluator("cell_volume")
      ->HasFieldChanged(S_next_.ptr(), name_);
  S_next_->GetFieldEvaluator("porosity")
      ->HasFieldChanged(S_next_.ptr(), name_);

  return false;
}


} // namespace
} // namespace
