/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

*/

#include "mpc_morphology_pk.hh"
#include "Mesh.hh"
#include "Mesh_MSTK.hh"


namespace Amanzi {

  Morphology_PK::Morphology_PK(Teuchos::ParameterList& pk_tree_or_fe_list,
                               const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                               const Teuchos::RCP<State>& S,
                               const Teuchos::RCP<TreeVector>& soln) :
    PK(pk_tree_or_fe_list, global_list, S, soln),
    PK_MPCSubcycled_ATS(pk_tree_or_fe_list, global_list, S, soln)
  {

     // Create verbosity object.
    vo_ = Teuchos::null;
    Teuchos::ParameterList vlist;
    vlist.sublist("verbose object") = plist_ -> sublist("verbose object");
    vo_ =  Teuchos::rcp(new VerboseObject("Morphology_PK", vlist)); 
    domain_ = plist_->get<std::string>("domain name", "domain");
    name_ = "morphology pk";

    Teuchos::Array<std::string> pk_order = plist_->get< Teuchos::Array<std::string> >("PKs order");
    
    
  }

// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double Morphology_PK::get_dt() {

  if (dt_MPC_ < 0) {
    double dt = Amanzi::PK_MPCSubcycled_ATS::get_dt();
    set_dt(dt);
    return dt;
  }else{
    return dt_MPC_;
  }
 
}


void Morphology_PK::Setup(const Teuchos::Ptr<State>& S){

  //passwd_ = "coupled_transport";  // owner's password
  //passwd_ = "state";  // owner's password


  dt_MPC_ = plist_->get<double>("dt MPC", 31557600);
  dt_sample_ = plist_->get<double>("dt sample", 259200);
  
  Amanzi::PK_MPCSubcycled_ATS::Setup(S);
  
  mesh_ = S->GetDeformableMesh(domain_);
  vertex_coord_key_ = Keys::getKey(domain_ , "vertex_coordinate");
  if (domain_ == "surface"){
    domain_3d_ = "surface_3d";
    domain_ss_ = "domain";
    vertex_coord_key_3d_ = Keys::getKey(domain_3d_ , "vertex_coordinate");
    vertex_coord_key_ss_ = Keys::getKey(domain_ss_ , "vertex_coordinate");
    mesh_3d_ = S->GetDeformableMesh(domain_ + "_3d");
    mesh_ss_ = S->GetDeformableMesh(domain_ss_);
  }
        
  // create storage for the vertex coordinates
  // we need to checkpoint those to be able to create
  // the deformed mesh after restart
  std::vector<AmanziMesh::Entity_kind> location(1);
  std::vector<int> num_dofs(1);
  std::vector<std::string> name(1);

  if (!S->HasField(vertex_coord_key_)){
    int dim = mesh_->space_dimension();
    location[0] = AmanziMesh::NODE;
    num_dofs[0] = dim;
    name[0] = "node";

    S->RequireField(vertex_coord_key_, "state")->SetMesh(mesh_)->SetGhosted()
      ->SetComponents(name, location, num_dofs);    
  }

  if (S->HasMesh(domain_3d_) && (!S->HasField(vertex_coord_key_3d_)) ){
    int dim = mesh_3d_->space_dimension();
    location[0] = AmanziMesh::NODE;
    num_dofs[0] = dim;
    name[0] = "node";

    S->RequireField(vertex_coord_key_3d_, "state")->SetMesh(mesh_3d_)->SetGhosted()
        ->SetComponents(name, location, num_dofs);        
  }

  if (S->HasMesh(domain_ss_) && (!S->HasField(vertex_coord_key_ss_)) ){
    int dim = mesh_ss_->space_dimension();
    location[0] = AmanziMesh::NODE;
    num_dofs[0] = dim;
    name[0] = "node";

    S->RequireField(vertex_coord_key_ss_, "state")->SetMesh(mesh_ss_)->SetGhosted()
        ->SetComponents(name, location, num_dofs);        
  }

  elevation_increase_key_ = Keys::getKey(domain_, "deformation");
  if (!S->HasField(elevation_increase_key_)){
    S->RequireField(elevation_increase_key_, "state")->SetMesh(mesh_)->SetGhosted(false)->SetComponent("cell", AmanziMesh::CELL, 1);
    Teuchos::ParameterList deform_plist;
    deform_plist.set("evaluator name", elevation_increase_key_);
    deform_eval_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(deform_plist));
    S->SetFieldEvaluator(elevation_increase_key_, deform_eval_);
  }

  int num_veg_species = plist_->get<int>("number of vegitation species", 1);

  
  Key biomass_key = Keys::getKey(domain_, "biomass");
  Key stem_density_key = Keys::getKey(domain_, "stem_density");
  Key stem_height_key = Keys::getKey(domain_, "stem_height");
  Key stem_diameter_key = Keys::getKey(domain_, "stem_diameter");
  Key plant_area_key = Keys::getKey(domain_, "plant_area");

  location[0] = AmanziMesh::CELL;
  num_dofs[0] = num_veg_species;
  name[0] = "cell";

  if (!S->HasField(biomass_key)){
    S->RequireField(biomass_key, biomass_key)->SetMesh(mesh_)->SetGhosted()
      ->SetComponents(name, location, num_dofs);
    S->RequireFieldEvaluator(biomass_key);
  }

  if (!S->HasField("msl"))
    S->RequireScalar("msl");
 
}

void Morphology_PK::Initialize(const Teuchos::Ptr<State>& S){

  Amanzi::PK_MPCSubcycled_ATS::Initialize(S);

  // initialize the vertex coordinate of existing meshes

  if (S->HasField(vertex_coord_key_))
    Initialize_MeshVertices_(S, mesh_, vertex_coord_key_);

  if (S->HasField(vertex_coord_key_3d_))
    Initialize_MeshVertices_(S, mesh_3d_, vertex_coord_key_3d_);

  if (S->HasField(vertex_coord_key_ss_))
    Initialize_MeshVertices_(S, mesh_ss_, vertex_coord_key_ss_);

  if (S->HasField(elevation_increase_key_)){
    S->GetFieldData(elevation_increase_key_, "state")->PutScalar(0.);
    S->GetField(elevation_increase_key_, "state")->set_initialized();
  }
    
  flow_pk_ = Teuchos::rcp_dynamic_cast<PK_BDF_Default> (sub_pks_[0]);
  sed_transport_pk_ = sub_pks_[1];


  const Epetra_MultiVector& dz = *S->GetFieldData(elevation_increase_key_)->ViewComponent("cell", false);

  dz_accumul_ = Teuchos::rcp(new Epetra_MultiVector(dz));
  dz_accumul_->PutScalar(0.);
}

void Morphology_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {
  sed_transport_pk_->CommitStep(t_old , t_new, S);

  S->set_time(t_new);
  // Key elev_key = Keys::readKey(*plist_, domain_, "elevation", "elevation");
  // Key slope_key = Keys::readKey(*plist_, domain_, "slope magnitude", "slope_magnitude");
    
  // bool chg = S_ -> GetFieldEvaluator(elev_key)->HasFieldChanged(S_.ptr(), elev_key);
  // if (chg){
  //   Teuchos::RCP<CompositeVector> elev =  S->GetFieldData(elev_key, elev_key);
  //   Teuchos::RCP<CompositeVector> slope = S->GetFieldData(slope_key, slope_key);
  //   Teuchos::RCP<CompositeVector> vc = S->GetFieldData(vertex_coord_key_ss_, "state");
  //   Teuchos::RCP<CompositeVector> dz = S->GetFieldData(elevation_increase_key_, "state");

  //   *elev  = *S_ -> GetFieldData(elev_key, elev_key);
  //   *slope = *S_ -> GetFieldData(slope_key, slope_key);
  //   *vc    = *S_ -> GetFieldData(vertex_coord_key_ss_, "state");
  // }

  // Key biomass_key = Keys::getKey(domain_, "biomass");
  // chg = S_ -> GetFieldEvaluator(biomass_key) -> HasFieldChanged(S_.ptr(), biomass_key);
  // if (chg)
    
  
}

// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool Morphology_PK::AdvanceStep(double t_old, double t_new, bool reinit) {
  bool fail = false;
  Teuchos::OSTab tab = vo_->getOSTab();
  
  double dt_step;
  if (dt_sample_ > 0) dt_step = std::min(S_->final_time() - S_->initial_time(), dt_sample_);
  else dt_step = S_->final_time() - S_->initial_time();

  
  Teuchos::RCP<Field_Scalar> msl_rcp = Teuchos::rcp_dynamic_cast<Field_Scalar>(S_->GetField("msl", "state"));
  msl_rcp->Compute(t_old);

  Key elev_key = Keys::readKey(*plist_, domain_, "elevation", "elevation");
  Epetra_MultiVector& dz = *S_next_->GetFieldData(elevation_increase_key_, "state")->ViewComponent("cell",false);
  dz.PutScalar(0.);

  flow_pk_ -> ResetTimeStepper(t_old);
  
  S_->set_intermediate_time(t_old);
  S_next_->set_intermediate_time(t_old);
  double dt_done = 0;
  double dt_next = flow_pk_ -> get_dt();
  double t_sample_end = t_old + dt_step;
  
  bool done = false;
  int ncycles = 0;

  while(!done){
    dt_next = flow_pk_ -> get_dt();   
    if (t_old + dt_done + dt_next > t_sample_end) {
      dt_next = t_sample_end - t_old - dt_done;
    }
    
    fail = true;
    while (fail){
      S_next_ -> set_time(t_old + dt_done + dt_next);
      S_ -> set_time(t_old + dt_done);
      S_inter_-> set_time(t_old + dt_done);
      fail = flow_pk_ -> AdvanceStep(t_old + dt_done, t_old + dt_done + dt_next, reinit);
      fail |= !flow_pk_->ValidStep();
      
      if (fail) {
        if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) *vo_->os()<<"Master step is failed\n";
        dt_next = flow_pk_ -> get_dt();
      }
      
    }

    master_dt_ = dt_next;
    //flow_pk_ -> CalculateDiagnostics(S_next_);
    flow_pk_ -> CommitStep(t_old  + dt_done, t_old + dt_done + dt_next, S_next_);

   
    slave_dt_ = sed_transport_pk_->get_dt(); 
    if (slave_dt_ > master_dt_) slave_dt_ = master_dt_;
    if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH)
      *vo_->os()<<"Slave dt="<<slave_dt_<<" Master dt="<<master_dt_<<"\n"; 
   
    fail = sed_transport_pk_->AdvanceStep(t_old + dt_done, t_old + dt_done + dt_next, reinit);

    if (fail){
      dt_next /= 2;
    }else{
      S_ -> set_intermediate_time(t_old + dt_done + dt_next);
      sed_transport_pk_->CommitStep(t_old + dt_done, t_old + dt_done + dt_next, S_next_);
      dt_done += dt_next;

      // we're done with this time step, copy the state
      *S_ = *S_next_;
      *S_inter_ = *S_next_;
      
    }
    ncycles ++;


    // check for done condition
    done = (std::abs(t_old + dt_done - t_sample_end) / (t_sample_end - t_old) < 0.1*min_dt_) || // finished the step
      (dt_next  < min_dt_); // failed
  }


  if ((dt_sample_ > 0)&&(dt_MPC_ >0)){
    double factor = dt_MPC_ / dt_sample_;   
    dz.Scale(factor);
  }


  dz_accumul_->Update(1, dz, 1);

  
  Update_MeshVertices_(S_next_.ptr() );

  
  bool chg = S_next_ -> GetFieldEvaluator(elev_key)->HasFieldChanged(S_next_.ptr(), elev_key);
  Epetra_MultiVector& elev_cell = *S_next_->GetFieldData(elev_key, elev_key)->ViewComponent("cell",false);
  for (int c=0; c<15; c++) std::cout<<c<<" "<<dz[0][c]<<" "<<(*dz_accumul_)[0][c]<<" "<<elev_cell[0][c] <<"\n";

  return fail;


  // // advance the slave, subcycling if needed
  // S_next_->set_intermediate_time(t_old);
  // bool done = false;
  // double dt_next = slave_dt_;
  // double dt_done = 0.;
  // int ncycles = 0;
  
  // while (!done) {
  //   // do not overstep
  //   if (t_old + dt_done + dt_next > t_new) {
  //     dt_next = t_new - t_old - dt_done;
  //   }
  //   // take the step
  //   fail = sed_transport_pk_->AdvanceStep(t_old + dt_done, t_old + dt_done + dt_next, reinit);
  //   ncycles ++;    
  //   if (fail) {
  //     // if fail, cut the step and try again
  //     dt_next /= 2;
  //   } else {
  //     // if success, commit the state and increment to next intermediate
  //     // -- etc: unclear if state should be commited or not?
  //     // set the intermediate time
  //     S_next_ -> set_intermediate_time(t_old + dt_done + dt_next);
  //     //S_next_ -> set_intermediate_time(t_old + dt_done + dt_next);
  //     sed_transport_pk_->CommitStep(t_old + dt_done, t_old + dt_done + dt_next, S_next_);
  //     //sed_transport_pk_->CommitStep(t_old + dt_done, t_old + dt_done + dt_next, S_next_);
  //     dt_done += dt_next;
  //   }
  //   // check for done condition
  //   done = (std::abs(t_old + dt_done - t_new) / (t_new - t_old) < 0.1*min_dt_) || // finished the step
  //       (dt_next  < min_dt_); // failed
  // }


  // if (std::abs(t_old + dt_done - t_new) / (t_new - t_old) < 0.1*min_dt_) {
  //   // done, success
  //   if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) *vo_->os()<<"Slave step is successful after "
  //                                                            <<ncycles <<" subcycles\n";
  //   fail = false;
  // } else {
  //   if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) *vo_->os()<<"Slave step is failed after "
  //                                                            <<ncycles <<" subcycles\n";
  //   fail = true;
  // }  
  
  return fail;

}

void Morphology_PK::Initialize_MeshVertices_(const Teuchos::Ptr<State>& S,
                                             Teuchos::RCP<const AmanziMesh::Mesh> mesh,
                                             Key vert_field_key){

  // spatial dimension
  int dim = mesh -> space_dimension();
  Amanzi::AmanziGeometry::Point coords(dim);
  // number of vertices
  int nV = mesh -> num_entities(Amanzi::AmanziMesh::NODE,
                                Amanzi::AmanziMesh::Parallel_type::OWNED);

  Epetra_MultiVector& vc = *S->GetFieldData(vert_field_key, "state")
      ->ViewComponent("node",false);

  // search the id of the mid point on the top
  for (int iV=0; iV<nV; iV++) {
    // get the coords of the node
    mesh -> node_get_coordinates(iV,&coords);
    for ( int s=0; s<dim; ++s ) {
      vc[s][iV] = coords[s];
    }
  }

  S -> GetFieldData(vert_field_key, "state") -> ScatterMasterToGhosted("node"); 
  S -> GetField(vert_field_key, "state") -> set_initialized();


}

void Morphology_PK::Update_MeshVertices_(const Teuchos::Ptr<State>& S){

  // spatial dimension
  int dim = mesh_ss_ -> space_dimension();
  Amanzi::AmanziGeometry::Point coords(dim);
  // number of vertices

  Epetra_MultiVector& vc = *S->GetFieldData(vertex_coord_key_ss_, "state")
    ->ViewComponent("node",true);
 
  
  const Epetra_MultiVector& dz = *S->GetFieldData(elevation_increase_key_)->ViewComponent("cell");

  int ncells = dz.MyLength();
  
  AmanziMesh::Entity_ID_List nodes, cells;
  double xyz[3];
  
  for (int c=0; c<ncells; c++) {

    AmanziMesh::Entity_ID domain_face;
    domain_face = mesh_->entity_get_parent(AmanziMesh::CELL, c);
    
    mesh_ss_ -> face_get_nodes(domain_face, &nodes);
    int nnodes = nodes.size();
    for (int i=0; i<nnodes; i++){
      mesh_ss_ -> node_get_coordinates(nodes[i], &coords);
      mesh_ss_ -> node_get_cells(nodes[i], Amanzi::AmanziMesh::Parallel_type::OWNED, &cells);
      int nsize = cells.size();
      double old = coords[2];

      coords[2] += dz[0][c] / nsize;
      vc[2][nodes[i]] += dz[0][c] / nsize;

      // coords[2] += 0.1 / nsize;
      // vc[2][nodes[i]] += 0.1 / nsize;

       // if (c==0){
       //   std::cout<<"COORDS "<<i<<" "<<old <<" "<<coords[2]<<" "<<vc[2][nodes[i]]<<"\n";
       // }
      
      mesh_ss_ -> node_set_coordinates(nodes[i], coords);
      
    }  
  }

  deform_eval_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(S -> GetFieldEvaluator(elevation_increase_key_)); 
  deform_eval_ -> SetFieldAsChanged(S.ptr());



  S -> GetFieldData(vertex_coord_key_, "state") -> ScatterMasterToGhosted("node"); 

}  


void Morphology_PK::FlowAnalyticalSolution_(const Teuchos::Ptr<State>& S, double time){

  Key mass_flux_key = "surface-mass_flux";
  Key depth_key = "surface-ponded_depth";
  
  Key owner;
  owner = S->GetField(mass_flux_key)->owner();
  Epetra_MultiVector& mass_flux_f = *S->GetFieldData(mass_flux_key, owner)->ViewComponent("face");
  owner = S->GetField(depth_key)->owner();
  Epetra_MultiVector& depth_c =  *S->GetFieldData(depth_key, owner)->ViewComponent("cell");

  Key elev_key = Keys::readKey(*plist_, domain_, "elevation", "elevation");
  const Epetra_MultiVector& elev_c = *S->GetFieldData(elev_key)->ViewComponent("cell");
  
  int ncells = depth_c.MyLength();
  int nfaces = mass_flux_f.MyLength();

  double AMP = 0.5;
  double tTide = 43200;
  double pi =  4 * std::atan(1.0);
  double wTide = 2*pi/tTide;
  double e_s = 0.3;
  double zmean = 0;

  for (int c=0; c<ncells; c++) zmean += elev_c[0][c];
  zmean *= 1./ncells;

  zmean = 0.3;
    
  double eta0  = -AMP*std::cos( wTide * time);

  double csi=-2.*(eta0 - zmean)/e_s;
  double psicsi, dsue_s, phismal;

  
  
  if (csi > 3 ) {
    psicsi=0.11053e-04;
    dsue_s=0.82777e-06;
    phismal=0.96894e-06;
  }else if(csi < -10) {
    psicsi=1.;
    dsue_s=(eta0 - zmean)/e_s;
    phismal=std::pow(dsue_s, 1.5);
  }else{
    double sign=1.;
    if(csi < 0) sign=-1.;
    double abscsi= std::abs(csi);
    double t=1./(1.+0.3275911*abscsi);
    double aux1=0.254829592*t - 0.284496736*(std::pow(t,2)) + 1.421413741*(std::pow(t,3));
    double aux2=-1.453152027*(std::pow(t,4))+1.061405429*(std::pow(t,5));
    double erfcsi=sign*(1.-(aux1+aux2)*exp(-abscsi*abscsi));

    psicsi=0.5*(1.-erfcsi);
    dsue_s=0.25*(exp(-csi*csi)/sqrt(pi)-csi*(1.-erfcsi));
    phismal=std::pow(dsue_s+0.27*sqrt(dsue_s)*exp(-2.*dsue_s), 1.5);
  }
  double depth=dsue_s*e_s;

  
  
  for (int c=0; c<ncells; c++){
    depth_c[0][c] = depth;   
  }

  // for (int f=0; f<nfaces; f++){
  //   AmanziGeometry::Point normal = S->GetMesh(domain_)->face_normal(f);
  //   if (std::abs(normal[1]) > 1e-4){
  //     AmanziGeometry::Point xf = S->GetMesh(domain_)->face_centroid(f);
  //     mass_flux_f[0][f] = -
  // }

  
  // Teuchos::RCP<FieldEvaluator> fm = S->GetFieldEvaluator(depth_key);
  
  
  

}


}  // namespace Amanzi
