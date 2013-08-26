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

  std::string mode_name = plist_.get<std::string>("deformation mode", "dVdt");
  if (mode_name == "dVdt") {
    mode_ = DEFORM_MODE_DVDT;
  } else if (mode_name == "thaw front") {
    mode_ = DEFORM_MODE_THAW_FRONT;
    deform_region_ = plist_.get<std::string>("deformation region");
    min_vol_frac_ = plist_.get<double>("minimum volume fraction");
  } else {
    Errors::Message mesg("Unknown deformation mode specified.");
    Exceptions::amanzi_throw(mesg);
  }
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
      S->RequireField("cell_volume_change", name_);
  cv_fac->Update(*rv_fac);

  switch(mode_) {
    case (DEFORM_MODE_DVDT): {
      // Create the deformation function
      Teuchos::ParameterList func_plist = plist_.sublist("deformation function");
      deform_func_ = Functions::CreateCompositeVectorFunction(func_plist, *cv_fac);
      break;
    }

    case (DEFORM_MODE_THAW_FRONT): {
      // Create storage for the initial face centroids
      S->RequireField("initial_face_height", name_)->SetMesh(mesh_)
          ->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 1);

      // Create storage for the initial cell volumes
      S->RequireField("initial_cell_volume", name_)->SetMesh(mesh_)
          ->SetComponent("cell", AmanziMesh::CELL, 1);

      // create the function to determine the front location
      Teuchos::ParameterList func_plist =
          plist_.sublist("thaw front function");
      FunctionFactory fac;
      thaw_front_func_ = Teuchos::rcp(fac.Create(func_plist));
      break;
    }

    default:
      ASSERT(0);
  }

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

  // create the operator
  def_matrix_ = Teuchos::rcp(new Operators::MatrixVolumetricDeformation(
      op_plist, mesh_));

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
  S->GetFieldData("cell_volume_change",name_)->PutScalar(0.);
  S->GetField("cell_volume_change",name_)->set_initialized();

  if (mode_ == DEFORM_MODE_THAW_FRONT) {
    {  // initialize the face centroid locations
      Epetra_MultiVector& face_height =
          *S->GetFieldData("initial_face_height",name_)
          ->ViewComponent("face",true);

      unsigned int nfaces_ghosted = face_height.MyLength();
      for (unsigned int f=0; f!=nfaces_ghosted; ++f) {
        face_height[0][f] = mesh_->face_centroid(f)[2];
      }

      S->GetField("initial_face_height",name_)->set_initialized();
    }

    { // initialize the cell volume
      Epetra_MultiVector& cv =
          *S->GetFieldData("initial_cell_volume",name_)
          ->ViewComponent("cell",false);

      unsigned int ncells = cv.MyLength();
      for (unsigned int c=0; c!=ncells; ++c) {
        cv[0][c] = mesh_->cell_volume(c);
      }

      S->GetField("initial_cell_volume",name_)->set_initialized();
    }

  }

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
      S_next_->GetFieldData("cell_volume_change", name_);

  Teuchos::RCP<AmanziMesh::Entity_ID_List> bottom_node_list =
      Teuchos::rcp(new AmanziMesh::Entity_ID_List());

  switch (mode_) {
    case (DEFORM_MODE_DVDT): {
      dcell_vol_vec->PutScalar(0.);
      deform_func_->Compute(T_mid, dcell_vol_vec.ptr());

      // scale by cell volume, this is fractional loss
      Epetra_MultiVector& dcell_vol_c = *dcell_vol_vec->ViewComponent("cell",false);
      for (int c=0; c!=dcell_vol_c.MyLength(); ++c) {
        dcell_vol_c[0][c] *= cv[0][c] * dT;
      }

      // collect a set of the fixed nodes
      std::vector<std::string> bottom_regions;
      Teuchos::ParameterList op_plist = plist_.sublist("global solve operator");
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
      for (std::set<AmanziMesh::Entity_ID>::const_iterator n=bottom_node_set.begin();
           n!=bottom_node_set.end(); ++n) {
        if (*n < nnodes_owned) bottom_node_list->push_back(*n);
      }

      break;
    }

    case (DEFORM_MODE_THAW_FRONT): {
      dcell_vol_vec->PutScalar(0.);

      const Epetra_MultiVector& face_height =
          *S_->GetFieldData("initial_face_height")->ViewComponent("face",true);
      const Epetra_MultiVector& cv0 =
          *S_->GetFieldData("initial_cell_volume")->ViewComponent("cell",false);
      const Epetra_MultiVector& cv_prev =
          *S_->GetFieldData("cell_volume")->ViewComponent("cell",false);

      double thaw_height = (*thaw_front_func_)(&ss);
      Epetra_MultiVector& dcell_vol_c = *dcell_vol_vec->ViewComponent("cell",false);

      double smallest_cv = 1.e10;
      double smallest_cv_frac = 1.;
      AmanziMesh::Entity_ID_List cells;
      mesh_->get_set_entities(deform_region_, AmanziMesh::CELL,
              AmanziMesh::OWNED, &cells);
      for (AmanziMesh::Entity_ID_List::const_iterator c=cells.begin();
           c!=cells.end(); ++c) {
        // determine up, down faces
        int my_up_n = -1;
        int my_down_n = -1;

        AmanziMesh::Entity_ID_List faces;
        std::vector<int> dirs;
        mesh_->cell_get_faces_and_dirs(*c, &faces, &dirs);
        double eps = 1.e-4;
        for (int n=0; n!=faces.size(); ++n) {
          if (mesh_->face_normal(faces[n],false,*c)[2] > eps) {
            ASSERT(my_up_n < 0);
            my_up_n = n;
          } else if (mesh_->face_normal(faces[n],false,*c)[2] < -eps) {
            ASSERT(my_down_n < 0);
            my_down_n = n;
          }
        }
        ASSERT(my_up_n >= 0);
        ASSERT(my_down_n >= 0);

        double cv_frac = std::min(std::max(
            (thaw_height - face_height[0][faces[my_down_n]]) /
            (face_height[0][faces[my_up_n]] - face_height[0][faces[my_down_n]]),
            min_vol_frac_), 1.);
        smallest_cv_frac = std::min(smallest_cv_frac, cv_frac);
        smallest_cv = std::min(smallest_cv, cv_frac * cv0[0][*c]);
        dcell_vol_c[0][*c] = cv_frac*cv0[0][*c] - cv_prev[0][*c];
      }

      std::cout << " Smallest cv_frac, cv = " << smallest_cv_frac << ", " << smallest_cv << std::endl;

      // set up the fixed list
      unsigned int nnodes = mesh_->num_entities(AmanziMesh::NODE,AmanziMesh::OWNED);
      for (unsigned int n=0; n!=nnodes; ++n) {
        AmanziGeometry::Point nc(3);
        mesh_->node_get_coordinates(n, &nc);
        if (nc[2] < thaw_height) bottom_node_list->push_back(n);
      }

      break;
    }


    default:
      ASSERT(0);
  }

  // calculate the nodal deformation
  Teuchos::RCP<CompositeVector> nodal_dz_vec =
      S_next_->GetFieldData("nodal_dz",name_);
  Teuchos::RCP<CompositeVector> rhs =
      Teuchos::rcp(new CompositeVector(*nodal_dz_vec));

  def_matrix_->Assemble(bottom_node_list.ptr());
  def_matrix_->InitializeInverse();
  def_matrix_->ApplyRHS(*dcell_vol_vec, rhs.ptr(), bottom_node_list.ptr());
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
    int nsurfnodes = surf_mesh_->num_entities(AmanziMesh::NODE,
            AmanziMesh::USED);

    Entity_ID_List surface_nodeids, surface3d_nodeids;
    AmanziGeometry::Point_List surface_newpos, surface3d_newpos;

    for (int i=0; i!=nsurfnodes; ++i) {
      // get the coords of the node
      AmanziMesh::Entity_ID pnode =
          surf_mesh_->entity_get_parent(AmanziMesh::NODE, i);
      AmanziGeometry::Point coord_domain(dim);
      mesh_->node_get_coordinates(pnode, &coord_domain);

      surface3d_nodeids.push_back(i);
      surface3d_newpos.push_back(coord_domain);

      // surface points are two dimensional
      AmanziGeometry::Point coord_surface(dim-1);
      for (int s=0; s!=dim-1; ++s) coord_surface[s] = coord_domain[s];
      surface_nodeids.push_back(i);
      surface_newpos.push_back(coord_surface);
    }
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
