/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Markus Berndt
           Daniil Svyatskiy

   Interface for the Volumetric Deformation PK.

   <ParameterList name="volumetric deformation">
   <Parameter name="PK model" type="string" value="Prescibed Mesh Deformation"/>
   <Parameter name="Deformation method" type="string" value="method name"/>
   </ParameterList>

   ------------------------------------------------------------------------- */
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "LinearOperatorFactory.hh"
#include "CompositeVectorFunctionFactory.hh"

#include "volumetric_deformation.hh"
#include "porosity_evaluator.hh"

namespace Amanzi {
namespace Deform {
     
using namespace Amanzi::AmanziMesh;

// RegisteredPKFactory<VolumetricDeformation> VolumetricDeformation::reg_("volumetric deformation");

  VolumetricDeformation::VolumetricDeformation(Teuchos::Ptr<State> S, const Teuchos::RCP<Teuchos::ParameterList>& plist,
        Teuchos::ParameterList& FElist,
        const Teuchos::RCP<TreeVector>& solution):
    PKDefaultBase(S, plist, FElist, solution),
    PKPhysicalBase(S, plist, FElist, solution) {
  poro_key_ = plist_->get<std::string>("porosity key","base_porosity");
  dt_ = plist_->get<double>("initial time step", 1);
  deform_value_ = 0.;

  // The deformation mode describes how to calculate new cell volume from a
  // provided function and the old cell volume.
  std::string mode_name = plist_->get<std::string>("deformation mode", "dVdt");
  if (mode_name == "dVdt") {
    deform_mode_ = DEFORM_MODE_DVDT;
  } else if (mode_name == "thaw front") {
    deform_mode_ = DEFORM_MODE_THAW_FRONT;
    deform_region_ = plist_->get<std::string>("deformation region");
    min_vol_frac_ = plist_->get<double>("minimum volume fraction");
  } else if (mode_name == "saturation") {
    deform_mode_ = DEFORM_MODE_SATURATION;
    deform_region_ = plist_->get<std::string>("deformation region");
    deform_value_ = plist_->get<double>("deformation porosity value", 0.3);
    //min_vol_frac_ = plist_->get<double>("minimum volume fraction");
    min_S_liq_ = plist_->get<double>("minimum liquid saturation", 0.3);
  } else {
    Errors::Message mesg("Unknown deformation mode specified.  Valid: [dVdt, thaw front].");
    Exceptions::amanzi_throw(mesg);
  }

  // The deformation strategy describes how to calculate nodal deformation
  // from cell volume change.
  std::string strategy_name = plist_->get<std::string>("deformation strategy",
          "global optimization");
  if (strategy_name == "global optimization") {
    strategy_ = DEFORM_STRATEGY_GLOBAL_OPTIMIZATION;
  } else if (strategy_name == "mstk implementation") {
    strategy_ = DEFORM_STRATEGY_MSTK;
  } else if (strategy_name == "average") {
    strategy_ = DEFORM_STRATEGY_AVERAGE;
  } else {
    Errors::Message mesg("Unknown deformation strategy specified. Valid: [global optimization, mstk implementation, average]");
    Exceptions::amanzi_throw(mesg);
  }

  // collect a set of the fixed nodes
  if (plist_->isParameter("bottom region")) {
    fixed_regions_.push_back(plist_->get<std::string>("bottom region"));
  } else {
    fixed_regions_ =
        plist_->get<Teuchos::Array<std::string> >("bottom regions").toVector();
  }

  fixed_region_type_ = plist_->get<std::string>("bottom region type", "node");
}

// -- Setup data
void VolumetricDeformation::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBase::setup(S);

  // save the meshes
  mesh_nc_ = S->GetDeformableMesh("domain");

  if (S->HasMesh("surface")) {
    surf_mesh_ = S->GetMesh("surface");
    surf3d_mesh_ = S->GetMesh("surface_3d");
    surf_mesh_nc_ = S->GetDeformableMesh("surface");
    surf3d_mesh_nc_ = S->GetDeformableMesh("surface_3d");
  }

  // create storage for primary variable, rock volume
  S->RequireField(key_, name_)->SetMesh(mesh_)->SetComponent("cell", AmanziMesh::CELL, 1);

  // Create storage and a function for cell volume change
  Teuchos::RCP<CompositeVectorSpace> cv_fac =  S->RequireField("cell_volume_change", name_);
  cv_fac->SetMesh(mesh_)->SetComponent("cell", AmanziMesh::CELL, 1);

  switch(deform_mode_) {
    case (DEFORM_MODE_DVDT): {
      // Create the deformation function
      Teuchos::ParameterList func_plist = plist_->sublist("deformation function");
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
      S->RequireField("integrated_cell_volume", name_)->SetMesh(mesh_)
          ->SetComponent("cell", AmanziMesh::CELL, 1);

      // create the function to determine the front location
      Teuchos::ParameterList func_plist =
          plist_->sublist("thaw front function");
      FunctionFactory fac;
      thaw_front_func_ = Teuchos::rcp(fac.Create(func_plist));
      break;
    }
    case (DEFORM_MODE_SATURATION): {
      // Create storage for the initial face centroids
      // S->RequireField("initial_face_height", name_)->SetMesh(mesh_)
      //     ->SetGhosted()->SetComponent("face", AmanziMesh::FACE, 1);

      // Create storage for the initial cell volumes
      // S->RequireField("initial_cell_volume", name_)->SetMesh(mesh_)
      //   ->SetComponent("cell", AmanziMesh::CELL, 1);
      // S->RequireField("integrated_cell_volume", name_)->SetMesh(mesh_)
      //   ->SetComponent("cell", AmanziMesh::CELL, 1);

      // create the function to determine the front location
      // Teuchos::ParameterList func_plist =
      //     plist_->sublist("thaw front function");
      // FunctionFactory fac;
      // thaw_front_func_ = Teuchos::rcp(fac.Create(func_plist));
      break;
    }

    default:
      ASSERT(0);
  }

  // create storage for the vertex coordinates
  // we need to checkpoint those to be able to create
  // the deformed mesh after restart
  int dim = mesh_->space_dimension();
  S->RequireField("vertex_coordinate_domain", name_)
      ->SetMesh(mesh_)->SetGhosted()
      ->SetComponent("node", AmanziMesh::NODE, dim);
  if (S->HasMesh("surface")) {
    S->RequireField("vertex_coordinate_surface_3d", name_)
        ->SetMesh(surf3d_mesh_)->SetGhosted()
        ->SetComponent("node", AmanziMesh::NODE, dim);
    S->RequireField("vertex_coordinate_surface", name_)
        ->SetMesh(surf_mesh_)->SetGhosted()
        ->SetComponent("node", AmanziMesh::NODE, dim-1);
  }

  S->RequireFieldEvaluator("cell_volume");
  S->RequireFieldEvaluator("porosity");
  S->RequireField("porosity")->SetMesh(mesh_)->SetGhosted()->AddComponent("cell",AmanziMesh::CELL,1);

  // Strategy-specific setup
  switch (strategy_) {
    case (DEFORM_STRATEGY_GLOBAL_OPTIMIZATION) : {
      // create the operator
      Teuchos::ParameterList op_plist = plist_->sublist("global solve operator");
      def_matrix_ = Teuchos::rcp(new Operators::MatrixVolumetricDeformation(
          op_plist, mesh_));

      if (op_plist.isSublist("Solver")) {
        Teuchos::ParameterList solver_plist = op_plist.sublist("Solver");
        AmanziSolvers::LinearOperatorFactory<CompositeMatrix,CompositeVector,
                            CompositeVectorSpace> fac;
        operator_ = fac.Create("deformation linear solver",solver_plist,def_matrix_);
      } else {
        operator_ = def_matrix_;
      }

      // create storage for the nodal deformation
      S->RequireField("nodal_dz", name_)->SetMesh(mesh_)->SetGhosted()
          ->SetComponent("node", AmanziMesh::NODE, 1);
      break;
    }
    case (DEFORM_STRATEGY_AVERAGE) : {
      // create storage for the nodal deformation, and count for averaging
      S->RequireField("nodal_dz", name_)->SetMesh(mesh_)->SetGhosted()
          ->SetComponent("node", AmanziMesh::NODE, 2);

      // create cell-based storage for deformation of the face above the cell
      S->RequireField("face_above_deformation", name_)->SetMesh(mesh_)
          ->SetGhosted()->SetComponent("cell", AmanziMesh::CELL, 1);
      break;
    }

    default: {}
  }
}


// -- Initialize owned (dependent) variables.
void VolumetricDeformation::initialize(const Teuchos::Ptr<State>& S) {
  PKPhysicalBase::initialize(S);

  //the PK's initial condition sets the initial porosity.  From this, we
  // calculate the actual initial condition, which is the rock volume.
  //  std::cout<<"name "<<name_<<" key_ "<<key_<<"\n";

  Epetra_MultiVector& base_poro = *S->GetFieldData(key_,name_)
      ->ViewComponent("cell",false);
  AmanziMesh::Entity_ID_List cells;
  mesh_->get_set_entities(deform_region_, AmanziMesh::CELL, AmanziMesh::OWNED, &cells);
  for (AmanziMesh::Entity_ID_List::const_iterator c=cells.begin(); c!=cells.end(); ++c) {
    base_poro[0][*c] = deform_value_;
  }
  

  // initialize the deformation
  S->GetFieldData("cell_volume_change",name_)->PutScalar(0.);
  S->GetField("cell_volume_change",name_)->set_initialized();

  switch (deform_mode_) {
    case (DEFORM_MODE_THAW_FRONT): {
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

      { 
        // initialize the cell volume
        Epetra_MultiVector& cv =
            *S->GetFieldData("initial_cell_volume",name_)
            ->ViewComponent("cell",false);

        unsigned int ncells = cv.MyLength();
        for (unsigned int c=0; c!=ncells; ++c) {
          cv[0][c] = mesh_->cell_volume(c);
        }
        *S->GetFieldData("integrated_cell_volume", name_)
            ->ViewComponent("cell",false) = cv;
        *S->GetFieldData("cell_volume", name_)
            ->ViewComponent("cell",false) = cv;

        S->GetField("integrated_cell_volume",name_)->set_initialized();
        S->GetField("initial_cell_volume",name_)->set_initialized();
        S->GetField("cell_volume",name_)->set_initialized();
      }
      break;
    }
  case (DEFORM_MODE_SATURATION): {
      { 
        // initialize the cell volume
        //S_next_->GetFieldEvaluator("cell_volume") -> HasFieldChanged(S_next_.ptr(), name_);
        Epetra_MultiVector& cv =
            *S->GetFieldData("cell_volume","cell_volume")
            ->ViewComponent("cell",false);

        unsigned int ncells = cv.MyLength();
        for (unsigned int c=0; c!=ncells; ++c) {
          cv[0][c] = mesh_->cell_volume(c);
        }

      }
      break;
    }
    default: {}
  }

  switch (strategy_) {
    case (DEFORM_STRATEGY_GLOBAL_OPTIMIZATION) : {
      // initialize the initial displacement to be zero
      S->GetFieldData("nodal_dz",name_)->PutScalar(0.);
      S->GetField("nodal_dz",name_)->set_initialized();
      break;
    }
    case (DEFORM_STRATEGY_AVERAGE) : {
      // initialize the initial displacement to be zero
      S->GetFieldData("nodal_dz",name_)->PutScalar(0.);
      S->GetField("nodal_dz",name_)->set_initialized();
      S->GetFieldData("face_above_deformation",name_)->PutScalar(0.);
      S->GetField("face_above_deformation",name_)->set_initialized();
      break;
    }
    default: {}
  }

  { // initialize the vertex coordinate to the current mesh
    int dim = mesh_->space_dimension();
    AmanziGeometry::Point coords(dim);
    int nnodes = mesh_->num_entities(Amanzi::AmanziMesh::NODE,
            Amanzi::AmanziMesh::OWNED);

    Epetra_MultiVector& vc = *S->GetFieldData("vertex_coordinate_domain",name_)
        ->ViewComponent("node",false);
    for (int iV=0; iV!=nnodes; ++iV) {
      // get the coords of the node
      mesh_->node_get_coordinates(iV,&coords);
      for (int s=0; s!=dim; ++s) vc[s][iV] = coords[s];
    }
    S->GetField("vertex_coordinate_domain",name_)->set_initialized();
  }

  if (S->HasMesh("surface")) {
    // initialize the vertex coordinates of the surface meshes
    int dim = surf_mesh_->space_dimension();
    AmanziGeometry::Point coords(dim);
    int nnodes = surf_mesh_->num_entities(Amanzi::AmanziMesh::NODE,
            Amanzi::AmanziMesh::OWNED);

    Epetra_MultiVector& vc = *S->GetFieldData("vertex_coordinate_surface",name_)
        ->ViewComponent("node",false);
    for (int iV=0; iV!=nnodes; ++iV) {
      // get the coords of the node
      surf_mesh_->node_get_coordinates(iV,&coords);
      for (int s=0; s!=dim; ++s) vc[s][iV] = coords[s];
    }
    S->GetField("vertex_coordinate_surface",name_)->set_initialized();
  }

  if (S->HasMesh("surface_3d")) {
    // initialize the vertex coordinates of the surface meshes
    int dim = surf3d_mesh_->space_dimension();
    AmanziGeometry::Point coords(dim);
    int nnodes = surf3d_mesh_->num_entities(Amanzi::AmanziMesh::NODE,
            Amanzi::AmanziMesh::OWNED);

    Epetra_MultiVector& vc = *S->GetFieldData("vertex_coordinate_surface_3d",name_)
        ->ViewComponent("node",false);
    for (int iV=0; iV!=nnodes; ++iV) {
      // get the coords of the node
      surf3d_mesh_->node_get_coordinates(iV,&coords);
      for (int s=0; s!=dim; ++s) vc[s][iV] = coords[s];
    }
    S->GetField("vertex_coordinate_surface_3d",name_)->set_initialized();
  }

}


bool VolumetricDeformation::advance(double dt) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_MEDIUM))
    *vo_->os() << "Advancing deformation PK from time " << S_->time() << " to "
               << S_next_->time() << " with step size " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;


  std::vector<double> ss(1,S_next_->time());
  double ss0 = S_->time();
  double dT = ss[0] - ss0;
  double T_mid = (ss[0] + ss0) / 2.;
  double min_height = 1e+23;
  double max_height = -1e+23;

  // Collect data from state
  Teuchos::RCP<CompositeVector> dcell_vol_vec =
      S_next_->GetFieldData("cell_volume_change", name_);
  dcell_vol_vec->PutScalar(0.);

  // Required only for DEFORM_STRATEGY_GLOBAL_OPTIMIZATION
  Teuchos::RCP<AmanziMesh::Entity_ID_List> fixed_node_list;

  // Calculate the change in cell volumes
  switch (deform_mode_) {
    case (DEFORM_MODE_SATURATION): {
      const Epetra_MultiVector& cv = *S_next_->GetFieldData("cell_volume")->ViewComponent("cell",true);
      const Epetra_MultiVector& s_liq =
        *S_->GetFieldData("saturation_liquid")->ViewComponent("cell",false);
      const Epetra_MultiVector& s_ice =
        *S_->GetFieldData("saturation_ice")->ViewComponent("cell",false);
      const Epetra_MultiVector& s_gas =
        *S_->GetFieldData("saturation_gas")->ViewComponent("cell",false);      
      const Epetra_MultiVector& temp =
          *S_->GetFieldData("temperature")->ViewComponent("cell",false);
      const Epetra_MultiVector& poro =
          *S_->GetFieldData("porosity")->ViewComponent("cell",false);
      const Epetra_MultiVector& base_poro =
          *S_->GetFieldData("base_porosity")->ViewComponent("cell",false);

      Epetra_MultiVector& dcell_vol_c = *dcell_vol_vec->ViewComponent("cell",false);
      int dim = mesh_->space_dimension();
      AmanziGeometry::Point centroid(dim);

      

      double min_porosity =  plist_->get<double>("minimum porosity", 0.5);
      double scl = plist_->get<double>("deformation scaling", 1.);
      double rel_comp_limit = plist_->get<double>("relative compresibility limit", 0.2);
        

      AmanziMesh::Entity_ID_List cells;
      mesh_->get_set_entities(deform_region_, AmanziMesh::CELL,
              AmanziMesh::OWNED, &cells);
      
      for (AmanziMesh::Entity_ID_List::const_iterator c=cells.begin(); c!=cells.end(); ++c) {

        centroid = mesh_->cell_centroid(*c);
        max_height = std::max(max_height, centroid[dim - 1]);

        double frac = 0.;


        if (s_liq[0][*c] > min_S_liq_ ){ // perform deformation if s_liq > min_S_liq_
          if ((poro[0][*c] - base_poro[0][*c])/base_poro[0][*c] < rel_comp_limit){ // perform deformation
                                                                            // if pressure have been relaxed enough
            frac = std::min((base_poro[0][*c] - min_porosity)/(1 - min_porosity), 
                             scl*(s_liq[0][*c] - min_S_liq_)*base_poro[0][*c]);
          }
        }
             
        dcell_vol_c[0][*c] = -frac*cv[0][*c];

        double soil_mass_vol = cv[0][*c]*(1 - base_poro[0][*c]);
        //         std::cout<<*c<<" "<<cv[0][*c]<<" "<<dcell_vol_c[0][*c]<<" frac "<<frac<<" poro "<<poro[0][*c]<<" ice "<<s_ice[0][*c]<<" liq "<<s_liq[0][*c]<<" gas "<<s_gas[0][*c]<< " soil vol "<<soil_mass_vol<<"\n";

      }

      if (strategy_ == DEFORM_STRATEGY_GLOBAL_OPTIMIZATION ||
          strategy_ == DEFORM_STRATEGY_MSTK) {
      //   // set up the fixed list
        fixed_node_list = Teuchos::rcp(new AmanziMesh::Entity_ID_List());
        AmanziMesh::Entity_ID_List nodes;
        mesh_->get_set_entities("bottom face", AmanziMesh::NODE,
                                AmanziMesh::OWNED, &nodes);
        for (AmanziMesh::Entity_ID_List::const_iterator n=nodes.begin();
             n!=nodes.end(); ++n) {                  
          fixed_node_list->push_back(*n);
        }
      }
      //exit(0);
      break;
    }
    default:
      ASSERT(0);
  }


  // Deform the subsurface mesh
  switch (strategy_) {
    case (DEFORM_STRATEGY_MSTK) : {
      // collect needed data, ghosted
      // -- cell vol
      Teuchos::RCP<const CompositeVector> cv_vec = S_->GetFieldData("cell_volume");
      const Epetra_MultiVector& cv = *cv_vec->ViewComponent("cell");
      // Teuchos::RCP<const CompositeVector> cv_vec_next = S_next_->GetFieldData("cell_volume");
      // const Epetra_MultiVector& cv_next = *cv_vec_next->ViewComponent("cell");



      Teuchos::RCP<const CompositeVector> poro_vec = S_->GetFieldData("porosity");
      const Epetra_MultiVector& poro = *poro_vec->ViewComponent("cell");

      const Epetra_MultiVector& s_liq =
        *S_->GetFieldData("saturation_liquid")->ViewComponent("cell",false);
      const Epetra_MultiVector& s_ice =
        *S_->GetFieldData("saturation_ice")->ViewComponent("cell",false);
      const Epetra_MultiVector& s_gas =
        *S_->GetFieldData("saturation_gas")->ViewComponent("cell",false);
    

      // -- dcell vol
      //dcell_vol_vec->ScatterMasterToGhosted("cell");
      const Epetra_MultiVector& dcell_vol_c =
          *dcell_vol_vec->ViewComponent("cell",true);

      // data needed in vectors
      int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::USED);
      int nnodes = mesh_->num_entities(AmanziMesh::NODE, AmanziMesh::USED);
      std::vector<double> target_cell_vols(ncells);
      std::vector<double> min_cell_vols(ncells);
      int dim = mesh_->space_dimension();
      AmanziGeometry::Point centroid(dim);
      min_height = max_height;

      for (int c=0; c!=ncells; ++c) {
        target_cell_vols[c] = cv[0][c] + dcell_vol_c[0][c];
        min_cell_vols[c] = (1 - poro[0][c])* cv[0][c] +  (poro[0][c]*s_ice[0][c]) * cv[0][c];        
        // min vol is rock vol + ice + a bit
        if (fabs(cv[0][c] - target_cell_vols[c])/cv[0][c] > 1e-4 ){
          //          std::cout<<c<<": "<<" "<<cv[0][c]<<" "<<target_cell_vols[c]<<" "<<min_cell_vols[c]<<"\n";
          centroid = mesh_->cell_centroid(c);
          min_height = std::min(min_height, centroid[dim - 1]);
          ASSERT(min_cell_vols[c] <= target_cell_vols[c]);
        }
        else {
          target_cell_vols[c] = -1.; //disregard these cells
        }

      }


      //std::cout<<"min_height "<<min_height<<"\n";
        Teuchos::RCP<AmanziMesh::Entity_ID_List>  below_node_list = Teuchos::rcp(new AmanziMesh::Entity_ID_List());
      for (unsigned int n=0; n!=nnodes; ++n) {
        AmanziGeometry::Point nc(3);
        mesh_->node_get_coordinates(n, &nc);
        if (nc[dim-1] < min_height) 
          below_node_list->push_back(n);
      }
      //exit(0);

      // call deform
      mesh_nc_->deform(target_cell_vols, min_cell_vols, *below_node_list, true);

      break;
    }
    default :
      ASSERT(0);
  }

  Teuchos::RCP<const CompositeVector> cv_vec = S_->GetFieldData("cell_volume");
  //cv_vec->ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& cv = *cv_vec->ViewComponent("cell",true);

  // now we have to adapt the surface mesh to the new volume mesh
  // extract the correct new coordinates for the surface from the domain
  // mesh and update the surface mesh accordingly
  if (surf_mesh_ != Teuchos::null) {
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
      int dim = mesh_->space_dimension();
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

  {  // update vertex coordinates in state (for checkpointing and error recovery)
    Epetra_MultiVector& vc =
        *S_next_->GetFieldData("vertex_coordinate_domain",name_)
        ->ViewComponent("node",false);
    int dim = mesh_->space_dimension();
    int nnodes = vc.MyLength();
    for (int i=0; i!=nnodes; ++i) {
      AmanziGeometry::Point coords(dim);
      mesh_->node_get_coordinates(i,&coords);
      for (int s=0; s!=dim; ++s) vc[s][i] = coords[s];
    }
  }

  if (S_next_->HasMesh("surface")) {
    // update vertex coordinates in state (for checkpointing and error recovery)
    Epetra_MultiVector& vc =
        *S_next_->GetFieldData("vertex_coordinate_surface",name_)
        ->ViewComponent("node",false);
    int dim = surf_mesh_->space_dimension();
    int nnodes = vc.MyLength();
    for (int i=0; i!=nnodes; ++i) {
      AmanziGeometry::Point coords(dim);
      surf_mesh_->node_get_coordinates(i,&coords);
      for (int s=0; s!=dim; ++s) vc[s][i] = coords[s];
    }
  }

  if (S_next_->HasMesh("surface_3d")) {
    // update vertex coordinates in state (for checkpointing and error recovery)
    Epetra_MultiVector& vc =
        *S_next_->GetFieldData("vertex_coordinate_surface_3d",name_)
        ->ViewComponent("node",false);
    int dim = surf3d_mesh_->space_dimension();
    int nnodes = vc.MyLength();
    for (int i=0; i!=nnodes; ++i) {
      AmanziGeometry::Point coords(dim);
      surf3d_mesh_->node_get_coordinates(i,&coords);
      for (int s=0; s!=dim; ++s) vc[s][i] = coords[s];
    }
  }

  // update cell volumes
  S_next_->GetFieldEvaluator("cell_volume") -> HasFieldChanged(S_next_.ptr(), name_);


  // update base porosity
  Teuchos::RCP<const CompositeVector> cv_vec_new = S_next_->GetFieldData("cell_volume");

  cv_vec->ScatterMasterToGhosted("cell");
  const Epetra_MultiVector& cv_new = *cv_vec_new->ViewComponent("cell",false);

  Teuchos::RCP<const CompositeVector> base_poro_vec_old = S_->GetFieldData(poro_key_);
  const Epetra_MultiVector& base_poro_old = *base_poro_vec_old->ViewComponent("cell");
  // Output
  Epetra_MultiVector& base_poro = *S_next_->GetFieldData(key_, name_)->ViewComponent("cell");

  int ncells = base_poro.MyLength();
  for (int c=0; c!=ncells; ++c) {
    base_poro[0][c] = 1 - (1. - base_poro_old[0][c]) * cv[0][c]/cv_new[0][c];
    // if (fabs(cv[0][c]/cv_new[0][c] - 1)>0.001){
    //   std::cout<<"volumes "<<c<<": "<<cv[0][c]<<" "<<cv_new[0][c]<<"\n";
    //   std::cout<<"porosity "<<base_poro[0][c] <<" "<<base_poro_old[0][c]<<"\n";
    // }
  }  

  S_next_->GetFieldEvaluator("porosity") -> HasFieldChanged(S_next_.ptr(), name_);

  // const Epetra_MultiVector& poros = *S_next_->GetFieldData("porosity")->ViewComponent("cell");

  // std::cout<<poros<<"\n";


  return false;
}


} // namespace
} // namespace
