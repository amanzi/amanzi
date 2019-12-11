/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
 * ATS
 *
 * License: see $ATS_DIR/COPYRIGHT
 * Author: Daniil Svyatsky, Xu Chonggang
 *
 * DOCUMENT ME:

 *
 *
 * ------------------------------------------------------------------------- */


#include "fates_pk.hh"
#include "bgc_simple.hh"
#include "vegetation.hh"

namespace Amanzi {
namespace BGC {


FATES_PK::FATES_PK(Teuchos::ParameterList& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& solution):
  PK(pk_tree, global_list,  S, solution),
  PK_Physical_Default(pk_tree, global_list,  S, solution),
  ncells_per_col_(-1)
{
  
  domain_surf_ = plist_->get<std::string>("surface domain name", "surface");

  Teuchos::ParameterList& FElist = S->FEList();

  //timestep size
  dt_ = plist_->get<double>("max time step", 1.e99);  
  dt_photosynthesis_ = plist_->get<double>("photosynthesis time step", 1800);
  dt_site_dym_ = plist_->get<double>("veg dynamics time step", 86400);
  surface_only_ = plist_->get<bool>("surface only", true);

}



  
void FATES_PK::Setup(const Teuchos::Ptr<State>& S){

  PK_Physical_Default::Setup(S);

  dt_ = plist_->get<double>("initial time step", 1.);

  // my mesh is the subsurface mesh, but we need the surface mesh, index by column, as well
  mesh_surf_ = S->GetMesh("surface");
  if (!surface_only_) mesh_=S->GetMesh();
  
  //  soil_part_name_ = plist_->get<std::string>("soil partition name");

  int iulog=10;
  int masterproc = 0;
  CFI_cdesc_t fatesdesc, clmdesc;
  int retval;

  if (!plist_->isParameter("fates parameter file")) {
    Errors::Message msg;
    msg << "No fates parameter file found in the parameter list for 'FATES'.\n";
    Exceptions::amanzi_throw(msg);
  }
  if (!plist_->isParameter("clm parameter file")) {
    Errors::Message msg;
    msg << "No clm parameter file found in the parameter list for 'FATES'.\n";
    Exceptions::amanzi_throw(msg);
  }
   
  std::string fates_file_str = plist_->get<std::string>("fates parameter file");
  char *fates_file_chr = &fates_file_str[0];
  std::string clm_file_str = plist_->get<std::string>("clm parameter file");
  char *clm_file_chr = &clm_file_str[0];

  retval = CFI_establish(&fatesdesc, fates_file_chr, CFI_attribute_other, CFI_type_char, fates_file_str.length(), 0, NULL);
  retval = CFI_establish(&clmdesc, clm_file_chr, CFI_attribute_other, CFI_type_char, clm_file_str.length(), 0, NULL); 
  fatessetmasterproc(&masterproc);
  fatessetinputfiles(&clmdesc, &fatesdesc);    
  fatesreadparameters();
  /*  ! Read in FATES parameter values early in the call sequence as well
    ! The PFT file, specifically, will dictate how many pfts are used
    ! in fates, and this will influence the amount of memory we
    ! request from the model, which is relevant in set_fates_global_elements()*/
  
  fatesreadpfts();

  /*   ------------------------------------------------------------------------
       Ask Fates to evaluate its own dimensioning needs.
       This determines the total amount of space it requires in its largest
       dimension.  We are currently calling that the "cohort" dimension, but
       it is really a utility dimension that captures the models largest
       size need.
       Sets:
       fates_maxElementsPerPatch
       fates_maxElementsPerSite (where a site is roughly equivalent to a column)
     
       (Note: fates_maxELementsPerSite is the critical variable used by CLM
       to allocate space)
     ------------------------------------------------------------------------*/

  set_fates_global_elements();
  
  // Get from FATES the total number of cohort size class bins output
  get_nlevsclass(&nlevsclass_);
  
  // requirements: primary variable
  S->RequireField(key_, name_)->SetMesh(mesh_surf_)
    ->SetComponent("cell", AmanziMesh::CELL, nlevsclass_);

  patchno_ = plist_->get<int>("number of patches", 1);
  nlevdecomp_ = plist_->get<int>("number of decomposition levels", 1);

  // met_decomp_key_ = Keys::getKey(domain_surf_,"decomp_cpools_met");
  // if (!S->HasField(met_decomp_key_)){
  //   S->RequireField(met_decomp_key_, name_)->SetMesh(mesh_surf_)
  //     ->SetComponent("cell", AmanziMesh::CELL, nlevdecomp_);
  // }

  // cel_decomp_key_ = Keys::getKey(domain_surf_,"decomp_cpools_cel");
  // if (!S->HasField(cel_decomp_key_)){
  //   S->RequireField(cel_decomp_key_, name_)->SetMesh(mesh_surf_)
  //     ->SetComponent("cell", AmanziMesh::CELL, nlevdecomp_);
  // }

  // lig_decomp_key_ = Keys::getKey(domain_surf_,"decomp_cpools_lig");
  // if (!S->HasField(lig_decomp_key_)){
  //   S->RequireField(lig_decomp_key_, name_)->SetMesh(mesh_surf_)
  //     ->SetComponent("cell", AmanziMesh::CELL, nlevdecomp_);
  // }
  
  precip_key_ = Keys::getKey(domain_surf_,"precipitation_rain");
  if (!S->HasField(precip_key_)){    
    S->RequireField(precip_key_, "state")->SetMesh(mesh_surf_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(precip_key_);
  }

  air_temp_key_ = Keys::getKey(domain_surf_,"air_temperature");
  if (!S->HasField(air_temp_key_)){    
    S->RequireField(air_temp_key_, "state")->SetMesh(mesh_surf_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(air_temp_key_);
  }

  humidity_key_ = Keys::getKey(domain_surf_,"relative_humidity");
  if (!S->HasField(humidity_key_)){    
    S->RequireField(humidity_key_, "state")->SetMesh(mesh_surf_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(humidity_key_);
  }

  wind_key_ = Keys::getKey(domain_surf_,"wind");
  if (!S->HasField(wind_key_)){    
    S->RequireField(wind_key_, "state")->SetMesh(mesh_surf_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(wind_key_);
  }

  co2a_key_ = Keys::getKey(domain_surf_,"co2a");
  if (!S->HasField(co2a_key_)){    
    S->RequireField(co2a_key_, "state")->SetMesh(mesh_surf_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(co2a_key_);
  }

  if (!surface_only_){
    poro_key_ = "base_porosity";
    if (!S->HasField(poro_key_)){    
      S->RequireField(poro_key_, "state")->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
      S->RequireFieldEvaluator(poro_key_);
    }
  }
  
}

  
void FATES_PK::Initialize(const Teuchos::Ptr<State>& S){

  PK_Physical_Default::Initialize(S);

  ncells_owned_ = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  site_.resize(ncells_owned_);

  t_photosynthesis_ = S->time();
  t_site_dym_ = S->time();
  
  // if (!S->GetField(met_decomp_key_, name_)->initialized()){
  //   S->GetFieldData(met_decomp_key_, name_)->PutScalar(0.0);
  //   S->GetField(met_decomp_key_, name_)->set_initialized();
  // }
  // if (!S->GetField(cel_decomp_key_, name_)->initialized()){
  //   S->GetFieldData(cel_decomp_key_, name_)->PutScalar(0.0);
  //   S->GetField(cel_decomp_key_, name_)->set_initialized();
  // }
  // if (!S->GetField(lig_decomp_key_, name_)->initialized()){
  //   S->GetFieldData(lig_decomp_key_, name_)->PutScalar(0.0);
  //   S->GetField(lig_decomp_key_, name_)->set_initialized();
  // }    


  
  if (surface_only_) ncells_per_col_ = 1;
  else{
    for (unsigned int col=0; col!=ncells_owned_; ++col) {
      int f = mesh_surf_->entity_get_parent(AmanziMesh::CELL, col);
      BGC::ColIterator col_iter(*mesh_, f);
      std::size_t ncol_cells = col_iter.size();
      if (ncells_per_col_ < 0) {
        ncells_per_col_ = ncol_cells;
      } else {
        AMANZI_ASSERT(ncol_cells == ncells_per_col_);
      }
    }
  }

  int array_size = ncells_per_col_*ncells_owned_;
  t_soil_.resize(array_size);
  sat_.resize(array_size);
  eff_poro_.resize(array_size);
  poro_.resize(array_size);
  suc_.resize(array_size);
  
  clump_ = 1;
  
  for (int i=0; i<ncells_owned_; ++i){
    site_[i].nlevbed = ncells_per_col_;
    site_[i].nlevdecomp = nlevdecomp_;
    site_[i].patchno = patchno_;
    site_[i].temp_veg24_patch = 273;
    site_[i].altmax_lastyear_indx_col = 1;
    site_[i].latdeg = plist_->get<double>("latitude");
    site_[i].londeg = plist_->get<double>("longitude");;
  }

  /* Preliminary initialization of FATES */
  init_ats_fates(&ncells_owned_, site_.data());

  double zi[6], z[5], dz[5], dzsoil_decomp[5];
  /* Define soil layers */
  zi[0] = 0;
  for (int i=0; i<ncells_per_col_; i++){
    zi[i+1] = zi[i] + 1;
    z[i] = 0.5*(zi[i+1] + zi[i]);
    dz[i] = 1;
  }  
  dzsoil_decomp[0] = 1;

  /* Initialize soil layers in FATES*/
  for (int i=0; i<ncells_owned_; i++){
    int s = i+1;    
    init_soil_depths(&clump_, &s,  &(site_[i]), zi, dz, z, dzsoil_decomp);
  }

  /* Init cold start of FATES */
  init_coldstart(&clump_);

  Epetra_MultiVector& biomass = *S->GetFieldData(key_, name_)->ViewComponent("cell", false);
  double* data_ptr;
  int data_dim;

  biomass.PutScalar(0.);
  biomass.ExtractView(&data_ptr, &data_dim);
  calculate_biomass(data_ptr, data_dim, nlevsclass_);
  S->GetField(key_, name_)->set_initialized();

}

double FATES_PK::get_dt(){

  dt_ = std::min(t_photosynthesis_ + dt_photosynthesis_ - S_inter_->time(), dt_photosynthesis_);
  dt_ = std::min(t_site_dym_ +  dt_site_dym_ - S_inter_->time(), std::min(dt_, dt_site_dym_));

  return dt_;

}


bool FATES_PK::AdvanceStep(double t_old, double t_new, bool reinit){


  double dtime = t_new - t_old;

  S_next_->GetFieldEvaluator(precip_key_)->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& precip_rain = *S_next_->GetFieldData(precip_key_)->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator(wind_key_)->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& wind = *S_next_->GetFieldData(wind_key_)->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator(humidity_key_)->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& humidity = *S_next_->GetFieldData(humidity_key_)->ViewComponent("cell", false);
  
  S_next_->GetFieldEvaluator(air_temp_key_)->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& air_temp = *S_next_->GetFieldData(air_temp_key_)->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator(co2a_key_)->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& co2a = *S_next_->GetFieldData(co2a_key_)->ViewComponent("cell", false);


  
  bool run_photo = false;
  bool run_veg_dym = false;

  if (fabs(t_new - (t_photosynthesis_ + dt_photosynthesis_)) < 1e-12*t_new) run_photo = true;
  if (fabs(t_new - (t_site_dym_ + dt_site_dym_)) <  1e-12*t_new) run_veg_dym = true;


  if (run_photo){

    for (unsigned int c=0; c<ncells_owned_; ++c){
      if (surface_only_){
        t_soil_[c] = air_temp[0][c];
        poro_[c] = 0.5;
        eff_poro_[c] = 0.5;
        sat_[c] = poro_[c];
        suc_[c] = -10;
      }else{

        S_next_->GetFieldEvaluator(poro_key_)->HasFieldChanged(S_next_.ptr(), name_);
        const Epetra_Vector& poro_vec = *(*S_next_->GetFieldData(poro_key_)->ViewComponent("cell", false))(0);        
        FieldToColumn_(c, poro_vec, poro_.data() + c*ncells_per_col_, ncells_per_col_);

      }

    }    

    int array_size = t_soil_.size();
    wrap_btran(&array_size, t_soil_.data(), poro_.data(), eff_poro_.data(), sat_.data(), suc_.data());

    PhotoSynthesisInput photo_input;
    // calculate fractional day length
    double t_days = S_inter_->time() / 86400.;
    int doy = std::floor(std::fmod(t_days, 365.25));
    if (doy == 0) doy = 365;

    
    const double& patm = *S_next_->GetScalarData("atmospheric_pressure");
    QSat qsat;


    
    for (unsigned int c=0; c<ncells_owned_; ++c){
      photo_input.dayl_factor = DayLength(site_[c].latdeg, doy);

      double es, esdT, qs, qsdT;
      qsat(air_temp[0][c], patm, &es, &esdT, &qs, &qsdT);
      photo_input.eair = es*humidity[0][c];

      double o2a = 209460.0;
      photo_input.cair = co2a[0][c] * patm * 1.e-6;
      photo_input.oair = o2a * patm * 1.e-6;
      photo_input.t_veg = air_temp[0][c];
      photo_input.tgcm = air_temp[0][c];
      photo_input.esat_tv = es;
      photo_input.eair = patm;      
    }

    t_photosynthesis_ = t_new;
    
  }

  
  std::vector<double> temp_veg24_patch(1);
  std::vector<double> prec24_patch(1);
  std::vector<double> rh24_patch(1);
  std::vector<double> wind24_patch(1);
  std::vector<double>  h2osoi_vol_col;

  if (run_veg_dym){
    for (int c=0; c<ncells_owned_; c++){
      int s=c+1;

      temp_veg24_patch[0] = air_temp[0][c];
      site_[c].temp_veg24_patch = air_temp[0][c];
      prec24_patch[0] = precip_rain[0][c];
      wind24_patch[0] = wind[0][c];
      rh24_patch[0] = humidity[0][c];


      if (surface_only_){
        h2osoi_vol_col.resize(1);
        h2osoi_vol_col[0] = 0.5;
      }else{

        h2osoi_vol_col.resize(ncells_per_col_);
        
        S_next_->GetFieldEvaluator(poro_key_)->HasFieldChanged(S_next_.ptr(), name_);
        const Epetra_Vector& poro_vec = *(*S_next_->GetFieldData(poro_key_)->ViewComponent("cell", false))(0);        
        FieldToColumn_(c, poro_vec, h2osoi_vol_col.data(), ncells_per_col_);

      }      
             
      dynamics_driv_per_site(&clump_, &s, &(site_[c]), &dtime,
                             h2osoi_vol_col.data(),
                             temp_veg24_patch.data(),
                             prec24_patch.data(),
                             rh24_patch.data(),
                             wind24_patch.data());
    
    }
    t_site_dym_ = t_new;
  }


  return false;
}

  

void FATES_PK::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S){

  Epetra_MultiVector& biomass = *S->GetFieldData(key_, name_)->ViewComponent("cell", false);
  double* data_ptr;
  int data_dim;

  biomass.PutScalar(0.);

  biomass.ExtractView(&data_ptr, &data_dim);

  calculate_biomass(data_ptr, data_dim, nlevsclass_);
  
}


// helper function for pushing field to column
void FATES_PK::FieldToColumn_(AmanziMesh::Entity_ID col, const Epetra_Vector& vec,
                               double* col_vec, int ncol) {
  ColIterator col_iter(*mesh_, mesh_surf_->entity_get_parent(AmanziMesh::CELL, col), ncol);
  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    col_vec[i] = vec[col_iter[i]];
  }
}

// helper function for collecting column dz and depth
void FATES_PK::ColDepthDz_(AmanziMesh::Entity_ID col,
                            Teuchos::Ptr<Epetra_SerialDenseVector> depth,
                            Teuchos::Ptr<Epetra_SerialDenseVector> dz) {
  AmanziMesh::Entity_ID f_above = mesh_surf_->entity_get_parent(AmanziMesh::CELL, col);
  ColIterator col_iter(*mesh_, f_above, ncells_per_col_);

  AmanziGeometry::Point surf_centroid = mesh_->face_centroid(f_above);
  AmanziGeometry::Point neg_z(3);
  neg_z.set(0.,0.,-1);

  for (std::size_t i=0; i!=col_iter.size(); ++i) {
    // depth centroid
    (*depth)[i] = surf_centroid[2] - mesh_->cell_centroid(col_iter[i])[2];

    // dz
    // -- find face_below
    AmanziMesh::Entity_ID_List faces;
    std::vector<int> dirs;
    mesh_->cell_get_faces_and_dirs(col_iter[i], &faces, &dirs);

    // -- mimics implementation of build_columns() in Mesh
    double mindp = 999.0;
    AmanziMesh::Entity_ID f_below = -1;
    for (std::size_t j=0; j!=faces.size(); ++j) {
      AmanziGeometry::Point normal = mesh_->face_normal(faces[j]);
      if (dirs[j] == -1) normal *= -1;
      normal /= AmanziGeometry::norm(normal);

      double dp = -normal * neg_z;
      if (dp < mindp) {
        mindp = dp;
        f_below = faces[j];
      }
    }

    // -- fill the val
    (*dz)[i] = mesh_->face_centroid(f_above)[2] - mesh_->face_centroid(f_below)[2];
    AMANZI_ASSERT( (*dz)[i] > 0. );
    f_above = f_below;
  }

}



  
}
}
