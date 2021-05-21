#include "Teuchos_XMLParameterListHelpers.hpp"

#include "mpc_surface_subsurface_helpers.hh"
#include "strong_mpc.hh"

#include "weak_mpc_semi_coupled_deform.hh"
#include "weak_mpc_semi_coupled_helper.hh"
#include <iomanip>


namespace Amanzi {

unsigned WeakMPCSemiCoupledDeform::flag_star = 0;
unsigned WeakMPCSemiCoupledDeform::flag_star_surf = 0;


WeakMPCSemiCoupledDeform::WeakMPCSemiCoupledDeform(Teuchos::ParameterList& pk_tree,
        const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
        const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& solution)
  : PK(pk_tree, global_plist, S, solution),
    MPC<PK>(pk_tree, global_plist, S, solution),
    sg_model_(false), dynamic_sg_model_(false)
{
  // grab the list of subpks
  auto subpks = plist_->get<Teuchos::Array<std::string> >("PKs order");
  std::string colname = subpks[subpks.size()-1];
  subpks.pop_back();

  KeyTriple col_triple;
  bool is_ds = Keys::splitDomainSet(colname, col_triple);
  if (!is_ds) {
    Errors::Message msg;
    msg << "WeakMPCSemiCoupledDeform subpk: \"" << colname << "\" should be a domain-set PK of the form column_*-NAME";
    Exceptions::amanzi_throw(msg);
  }

  // add for the various columns based on GIDs of the surface system
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = S->GetMesh("surface");
  int ncols = surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int i=0; i!=ncols; ++i) {
    int gid = surf_mesh->cell_map(false).GID(i);
    std::stringstream domain_name_stream;
    domain_name_stream << std::get<0>(col_triple) << "_" << gid;
    subpks.push_back(Keys::getKey(domain_name_stream.str(), std::get<2>(col_triple)));
  }

  numPKs_ = subpks.size();

  PKFactory pk_factory;
  // -- create the star pk
  {
    // -- create the solution vector
    Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector());
    solution_->PushBack(pk_soln);
    // -- create the PK
    Teuchos::RCP<PK> pk = pk_factory.CreatePK(subpks[0], pk_tree, global_list_, S, pk_soln);
    sub_pks_.push_back(pk);
  }

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // -- create the lifted PKs
  int npks = subpks.size();
  for (int lcv_i=1; lcv_i!=npks; ++lcv_i) {
    // -- create the solution vector
    Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector());
    solution_->PushBack(pk_soln);

    // -- create the PK
    Teuchos::RCP<PK> pk = pk_factory.CreatePK(subpks[lcv_i], pk_tree, global_list_, S, pk_soln);
    sub_pks_.push_back(pk);

    // -- check IC -- restart columns on their respective processors
    Teuchos::Array<std::string> pkorder1 = global_list_->sublist("PKs").sublist(subpks[lcv_i]).get<Teuchos::Array<std::string> >("PKs order");
    Teuchos::Array<std::string> pkorder2 = global_list_->sublist("PKs").sublist(pkorder1[1]).get<Teuchos::Array<std::string> >("PKs order");

    for (int lcv_j = 0; lcv_j < pkorder2.size(); lcv_j++) {
      Teuchos::ParameterList& sub_list = global_list_->sublist("PKs").sublist(pkorder1[0]).sublist("initial condition");
      if (sub_list.isParameter("restart files, cycles")) {
        Teuchos::Array<std::string> restart = sub_list.get<Teuchos::Array<std::string> >("restart files, cycles");
        std::stringstream res_file;

        if (restart[0].rfind("/") != restart[0].size()-1) restart[0] += "/";
        res_file << restart[0] << "checkpoint_" << rank << "_" << restart[1] << ".h5";
        sub_list.set("restart file", res_file.str());
      }

      if (global_list_->sublist("PKs").sublist(pkorder2[lcv_j]).isSublist("initial condition")) {
        Teuchos::ParameterList& seb_list = global_list_->sublist("PKs").sublist(pkorder2[lcv_j]).sublist("initial condition");

        if (seb_list.isParameter("restart files, cycles")) {
          Teuchos::Array<std::string> restart = seb_list.get<Teuchos::Array<std::string> >("restart files, cycles");
          std::stringstream res_file;
          if (restart[0].rfind("/") != restart[0].size()-1) restart[0] += "/";
          res_file << restart[0] << "checkpoint_" << rank << "_" << restart[1] << ".h5";
          seb_list.set("restart file", res_file.str());
        }

      } else if (global_list_->sublist("PKs").sublist(pkorder2[lcv_j]).isParameter("PKs order")) {
        Teuchos::Array<std::string> pkorder3 = global_list_->sublist("PKs").sublist(pkorder2[lcv_j]).get<Teuchos::Array<std::string> >("PKs order");

        for (int lcv_k=0; lcv_k < pkorder3.size(); lcv_k++) {
          Teuchos::ParameterList& seb_list = global_list_->sublist("PKs").sublist(pkorder3[lcv_k]).sublist("initial condition");
          if (seb_list.isParameter("restart files, cycles")) {
            Teuchos::Array<std::string> restart = seb_list.get<Teuchos::Array<std::string> >("restart files, cycles");

            std::stringstream res_file;
            if (restart[0].rfind("/") != restart[0].size()-1) restart[0] += "/";
            res_file << restart[0] << "checkpoint_" << rank << "_" << restart[1] << ".h5";
            seb_list.set("restart file", res_file.str());
          }
        }
      } //PSS end
    }
  }
}



// must communicate dts since columns are serial
double WeakMPCSemiCoupledDeform::get_dt()
{
  double dt = 1.0e99;
  for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    dt = std::min<double>(dt,(*pk)->get_dt());
  }

  double dt_local = dt;
  S_->GetMesh("surface")->get_comm()->MinAll(&dt_local, &dt, 1);
  return dt;
}

// -----------------------------------------------------------------------------
// Set timestep for sub PKs
// -----------------------------------------------------------------------------
void WeakMPCSemiCoupledDeform::set_dt( double dt)
{
  for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->set_dt(dt);
  }
};



// -----------------------------------------------------------------------------
// Set up each PK
// -----------------------------------------------------------------------------
void
WeakMPCSemiCoupledDeform::Setup(const Teuchos::Ptr<State>& S)
{
  S->AliasMesh("surface", "surface_star");
  MPC<PK>::Setup(S);

  coupling_key_ = plist_->get<std::string>("coupling key"," ");
  subcycle_key_ = plist_->get<bool>("subcycle",false);

  // by default sg_model_ is false
  if (S->FEList().isSublist("surface_star-depression_depth")){
    sg_model_ = true;
    if (S->FEList().isSublist("surface_star-thaw_depth")){
      dynamic_sg_model_= true; }
    else{
      dynamic_sg_model_ = false; }
  } else {
    sg_model_ = false;
  }
  AMANZI_ASSERT(!(coupling_key_.empty()));
};


// surface_column cells are initialized from the subsurface columns, and then surface_star is initialized from surface cells, so order is important
void
WeakMPCSemiCoupledDeform::Initialize(const Teuchos::Ptr<State>& S)
{
  MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
  ++pk;
  for (; pk!=sub_pks_.end(); ++pk){
    (*pk)->Initialize(S);
  }
  MPC<PK>::Initialize(S);
}

bool
WeakMPCSemiCoupledDeform::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = false;
  if (coupling_key_ == "surface subsurface system: columns"){
    fail = CoupledSurfSubsurfColumns(t_old, t_new, reinit);
  }
  else if(coupling_key_ == "surface subsurface system: 3D"){
    // can be used for 3D subsurface and surface_star coupling as well.
  }
  return fail;
};

bool
WeakMPCSemiCoupledDeform::CoupledSurfSubsurfColumns(double t_old, double t_new, bool reinit)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  bool fail = false;
  double M_ = 0.0180153; //molar mass

  //ensure the star solution is marked as changed when the subsurface columns fail
  if (flag_star || flag_star_surf) {
    flag_star = 0;
    flag_star_surf=0;

    Teuchos::RCP<PK_BDF_Default> pk_sfstar = Teuchos::rcp_dynamic_cast<PK_BDF_Default>(sub_pks_[0]);
    AMANZI_ASSERT(pk_sfstar.get());
    pk_sfstar->ChangedSolution();
  }

  fail = sub_pks_[0]->AdvanceStep(t_old, t_new, reinit);
  int nfailed_surf = 0;
  if (fail)
    nfailed_surf++;

  int nfailed_local_sf = nfailed_surf;
  S_->GetMesh("surface")->get_comm()->SumAll(&nfailed_local_sf, &nfailed_surf, 1);

  if (vo_->os_OK(Teuchos::VERB_LOW))
    *vo_->os() << "surface_star system failed: "<< ((nfailed_surf >0) ? 1 : 0) << std::endl;

  if (nfailed_surf > 0) {
    flag_star_surf=1;
    return true;
  }

  //copying surface_star (2D) data (pressures/temperatures) to column surface (1D-cells)[all the surf column cells get updates]
  const Epetra_MultiVector& surfstar_temp = *S_next_->GetFieldData("surface_star-temperature")
    ->ViewComponent("cell", false);

  unsigned int size_temp = surfstar_temp.MyLength();
  AMANZI_ASSERT(size_temp == numPKs_ -1); // check if the subsurface columns are equal to the surface cells

  //copying pressure
  if (!sg_model_) {
    const Epetra_MultiVector& surfstar_pres = *S_next_->GetFieldData("surface_star-pressure")->ViewComponent("cell", false);
    for (unsigned c = 0; c < size_temp; c++) {
      if (surfstar_pres[0][c] > 101325.00) {
        std::stringstream name;
        int id = S_->GetMesh("surface")->cell_map(false).GID(c);
        name << "surface_column_" << id;
        Epetra_MultiVector& surf_pres = *S_inter_->GetFieldData(Keys::getKey(name.str(),"pressure"),
                S_inter_->GetField(Keys::getKey(name.str(),"pressure"))->owner())->ViewComponent("cell", false);
        surf_pres[0][0] = surfstar_pres[0][c];
      }
    }

  } else {
    const Epetra_MultiVector& vol_pd = *S_next_->GetFieldData("surface_star-volumetric_ponded_depth")
      ->ViewComponent("cell", false);
    const Epetra_MultiVector& molar_dens_liq = *S_next_->GetFieldData("surface_star-mass_density_liquid")
      ->ViewComponent("cell", false);
    const Epetra_Vector& gravity = *S_->GetConstantVectorData("gravity");
    double gz = -gravity[2];

    for (unsigned c=0; c<size_temp; c++) {
      double pres = vol_pd[0][c] * molar_dens_liq[0][c] * gz + 101325.0; // convert volumetric head to pressure

      if (pres > 101325.0) {
        std::stringstream name;
        int id = S_->GetMesh("surface")->cell_map(false).GID(c);

        name << "surface_column_" << id;
        Epetra_MultiVector& surf_pres = *S_inter_->GetFieldData(Keys::getKey(name.str(),"pressure"),
								S_inter_->GetField(Keys::getKey(name.str(),"pressure"))->owner())
          ->ViewComponent("cell", false);
        surf_pres[0][0] = pres;
      }
    }
  }

  //copying temperatures
  for (unsigned c = 0; c < size_temp; c++) {
    std::stringstream name, name_ss;
    int id = S_->GetMesh("surface")->cell_map(false).GID(c);
    name << "surface_column_" << id;
    name_ss << "column_" << id;

    Epetra_MultiVector& surf_temp = *S_inter_->GetFieldData(Keys::getKey(name.str(),"temperature"),
            S_inter_->GetField(Keys::getKey(name.str(),"temperature"))->owner())->ViewComponent("cell", false);

    surf_temp[0][0] = surfstar_temp[0][c];

    CopySurfaceToSubsurface(*S_inter_->GetFieldData(Keys::getKey(name.str(),"pressure")),
                            S_inter_->GetFieldData(Keys::getKey(name_ss.str(),"pressure"),
                                    S_inter_->GetField(Keys::getKey(name_ss.str(),"pressure"))->owner()).ptr());
    CopySurfaceToSubsurface(*S_inter_->GetFieldData(Keys::getKey(name.str(),"temperature")),
                            S_inter_->GetFieldData(Keys::getKey(name_ss.str(),"temperature"),
                                    S_inter_->GetField(Keys::getKey(name_ss.str(),"temperature"))->owner()).ptr());
  }

  // successful surface_star system updates columns -- changes previous solution.
  for (int i=1; i<numPKs_; i++) {
    Teuchos::RCP<PK> pk_domain = Teuchos::rcp_dynamic_cast<PK>(sub_pks_[i]);
    AMANZI_ASSERT(pk_domain.get());
    pk_domain->ChangedSolutionPK(S_inter_.ptr());
  }

  // FIX ME, need a comm here! --etc
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int nfailed = 0;
  double t0 = S_inter_->time();
  double t1 = S_next_->time();

  // skip the first
  auto sub_pk = sub_pks_.begin();
  ++sub_pk;

  // Loop over all subsurface PKs (columns)
  for (auto pk = sub_pk; pk!=sub_pks_.end(); ++pk){
    std::stringstream name_ss;
    int c = pk - sub_pks_.begin() - 1;
    name_ss << "column_" << S_->GetMesh("surface")->cell_map(false).GID(c);

    if (!subcycle_key_) {
      bool c_fail = (*pk)->AdvanceStep(t_old, t_new, reinit);
      if (vo_->os_OK(Teuchos::VERB_LOW))
        *vo_->os() << name_ss.str() << "failed? "<< ((c_fail) ? 1 : 0) <<" rank: "<<rank<<std::endl ;
      if (c_fail) nfailed++;

    } else {
      //subcycling needs to be fixed and cleaned.
      std::stringstream name;
      int id = S_->GetMesh("surface")->cell_map(false).GID(c);
      name << "surface_column_" << id;

      double loc_dt = 0; //revisit dt;
      bool done = false;
      double t = 0;
      bool cyc_flag  = false;
      // S_inter_->set_time(t0+t);
      // S_next_->set_time(t0 + t + loc_dt);
      // lets put dt=0 for a while---revisit
      double dt =0;
      while (!done) {
        bool fail_pk = false; //----revisit= (*pk)->advance(loc_dt);
        //fail_pk |= !(*pk)->valid_step();

        if (!fail_pk) {
          //--revisit (*pk)->commit_state(loc_dt, S_next_);
          t = t + loc_dt;
          double loc_dt_old = loc_dt;

          double loc_dt_new = (*pk)->get_dt();
          if (dt - (t+loc_dt_new) < 1.)
            loc_dt_new = dt -t;
          else if ( t + loc_dt_new < dt && (t+2*loc_dt_new > dt ) )
            loc_dt_new = 0.5*(dt - t);

          loc_dt = std::min<double>(loc_dt_new, dt - t);

          if (loc_dt <= 1.e-10) {
            done = true;
            if (cyc_flag) S_inter_->set_time(t0);

          } else {
            UpdateIntermediateStateParameters(S_next_, S_inter_,id);

            Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe1 =
              Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>
              (S_inter_->GetFieldEvaluator(Keys::getKey(name.str(),"water_source_temperature")));

            pfe1->SetFieldAsChanged(S_inter_.ptr());

            Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe2 =
              Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>
              (S_inter_->GetFieldEvaluator(Keys::getKey(name.str(),"conducted_energy_source")));

            pfe2->SetFieldAsChanged(S_inter_.ptr());

            Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe3 =
              Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>
              (S_inter_->GetFieldEvaluator(Keys::getKey(name.str(),"water_source")));

            pfe3->SetFieldAsChanged(S_inter_.ptr());

            Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe4 =
              Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>
              (S_inter_->GetFieldEvaluator(Keys::getKey(name_ss.str(),"water_source")));

            pfe4->SetFieldAsChanged(S_inter_.ptr());

            Teuchos::RCP<PK_PhysicalBDF_Default> pk_domain =
              Teuchos::rcp_dynamic_cast<PK_PhysicalBDF_Default>(sub_pks_[c+1]);
            AMANZI_ASSERT(pk_domain.get()); // make sure the pk_domain is not empty
            pk_domain->ChangedSolution(S_inter_.ptr());

            S_inter_->set_time(t0+t);
            S_next_->set_time( t0 + t + loc_dt);
          }

        } else {
          cyc_flag = true;

          UpdateNextStateParameters(S_next_, S_inter_, id);

          Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe1 =
            Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>
            (S_next_->GetFieldEvaluator(Keys::getKey(name.str(),"water_source_temperature")));

          pfe1->SetFieldAsChanged(S_next_.ptr());

          Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe2 =
            Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>
            (S_next_->GetFieldEvaluator(Keys::getKey(name.str(),"conducted_energy_source")));

          pfe2->SetFieldAsChanged(S_next_.ptr());

          Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe3 =
            Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>
            (S_next_->GetFieldEvaluator(Keys::getKey(name.str(),"water_source")));

          pfe3->SetFieldAsChanged(S_next_.ptr());

          Teuchos::RCP<PrimaryVariableFieldEvaluator> pfe4 =
            Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>
            (S_next_->GetFieldEvaluator(Keys::getKey(name_ss.str(),"water_source")));

          pfe4->SetFieldAsChanged(S_next_.ptr());

          Teuchos::RCP<PK_PhysicalBDF_Default> pk_domain =
            Teuchos::rcp_dynamic_cast<PK_PhysicalBDF_Default>(sub_pks_[c+1]);
          AMANZI_ASSERT(pk_domain.get()); // make sure the pk_domain is not empty
          pk_domain->ChangedSolution(S_next_.ptr());

          loc_dt = (*pk)->get_dt();
          S_inter_->set_time(t0+t);
          // revisit S_next_->set_time(t0 + t + loc_dt);
        }

      }
    }
  }

  // FIX ME -- use a real comm! --etc
  MPI_Barrier(MPI_COMM_WORLD);

  int nfailed_local = nfailed;
  S_->GetMesh("surface")->get_comm()->SumAll(&nfailed_local, &nfailed, 1);

  if (nfailed == 0){
    Epetra_MultiVector& surfstar_p = *S_next_->GetFieldData("surface_star-pressure",
            S_inter_->GetField("surface_star-pressure")->owner())->ViewComponent("cell", false);
    Epetra_MultiVector& surfstar_t = *S_next_->GetFieldData("surface_star-temperature",
            S_inter_->GetField("surface_star-temperature")->owner())->ViewComponent("cell", false);
    Epetra_MultiVector& surfstar_wc = *S_next_->GetFieldData("surface_star-water_content",
            S_inter_->GetField("surface_star-water_content")->owner())->ViewComponent("cell", false);


    if (!sg_model_) {
      for (unsigned c=0; c<size_temp; c++) {
        std::stringstream name;
        int id = S_->GetMesh("surface")->cell_map(false).GID(c);

        name << "surface_column_" << id;
        const Epetra_MultiVector& surf_p = *S_next_->GetFieldData(Keys::getKey(name.str(),"pressure"))
          ->ViewComponent("cell", false);
        const Epetra_MultiVector& surf_wc = *S_next_->GetFieldData(Keys::getKey(name.str(),"water_content"))
          ->ViewComponent("cell", false);
        if (surf_p[0][0] > 101325.0) {
          surfstar_p[0][c] = surf_p[0][0];
          surfstar_wc[0][c] = surf_wc[0][0];
        } else {
          surfstar_p[0][c] = 101325.00;
        }
      }

    } else {
      const Epetra_MultiVector& delta_max_v = *S_next_->GetFieldData("surface_star-maximum_ponded_depth")
        ->ViewComponent("cell", false);
      const Epetra_MultiVector& delta_ex_v = *S_next_->GetFieldData("surface_star-excluded_volume")
        ->ViewComponent("cell", false);

      const Epetra_Vector& gravity = *S_->GetConstantVectorData("gravity");
      double gz = -gravity[2];
      const double& p_atm = *S_->GetScalarData("atmospheric_pressure");

      for (unsigned c=0; c<size_temp; c++) {
        std::stringstream name;
        int id = S_->GetMesh("surface")->cell_map(false).GID(c);
        name << "surface_column_" << id;
        const Epetra_MultiVector& pd = *S_next_->GetFieldData(Keys::getKey(name.str(),"ponded_depth"))
          ->ViewComponent("cell", false);
        const Epetra_MultiVector& surf_wc = *S_next_->GetFieldData(Keys::getKey(name.str(),"water_content"))
          ->ViewComponent("cell", false);

        const Epetra_MultiVector& cv = *S_next_->GetFieldData(Keys::getKey(name.str(),"cell_volume"))
          ->ViewComponent("cell", false);

        const Epetra_MultiVector& molar_dens_liq = *S_next_->GetFieldData(Keys::getKey(name.str(),"mass_density_liquid"))
          ->ViewComponent("cell", false);

        if (pd[0][0] > 0) {
          double delta = FindVolumetricHead(pd[0][0], delta_max_v[0][c],delta_ex_v[0][c]);
          double pres = delta * molar_dens_liq[0][0] * gz + p_atm;
          surfstar_p[0][c] = pres;

          double vpd = 0;
          if (delta <= delta_max_v[0][c]) {
            vpd = std::pow(delta,2) * (2*delta_max_v[0][c] - 3*delta_ex_v[0][c]) / std::pow(delta_max_v[0][c],2)
              + std::pow(delta,3) * (2*delta_ex_v[0][c] - delta_max_v[0][c]) / std::pow(delta_max_v[0][c],3);
            // later call get the volumetric head given ponded depth to fix this
          } else {
            vpd = delta - delta_ex_v[0][c];
          }

          double vpd_pres = vpd *molar_dens_liq[0][0] * gz + p_atm;
          surfstar_wc[0][c] = (vpd_pres - p_atm)/ (gz * M_);
          surfstar_wc[0][c] *= cv[0][0];
        } else {
          surfstar_p[0][c] = 101325.0;
        }
      }
    }

    for (unsigned c=0; c<size_temp; c++) {
      std::stringstream name;
      int id = S_->GetMesh("surface")->cell_map(false).GID(c);
      name << "surface_column_" << id;
      const Epetra_MultiVector& surf_t = *S_next_->GetFieldData(Keys::getKey(name.str(),"temperature"))->ViewComponent("cell", false);
      surfstar_t[0][c] = surf_t[0][0];
    }

    // Mark surface_star-pressure evaluator as changed.
    Teuchos::RCP<PK_BDF_Default> pk_surf =
      Teuchos::rcp_dynamic_cast<PK_BDF_Default>(sub_pks_[0]);
    AMANZI_ASSERT(pk_surf.get());
    pk_surf->ChangedSolution();

    MPC<PK>::SubPKList::iterator pk1 = sub_pks_.begin();

    if (subcycle_key_) (*pk1)->CommitStep(t_old, t_new, S_next_);
  }

  if (nfailed > 0) {
    flag_star = 1;
    return true;
  } else {
    flag_star = 1;
    return false;
  }
}


double
WeakMPCSemiCoupledDeform::FindVolumetricHead(double d, double delta_max, double delta_ex)
{
  // refactor this to use a proper root finding algorithm, assuming that is what is begin done here? --etc
  double a = (2*delta_ex - delta_max) / std::pow(delta_max,3);
  double b = (2*delta_max - 3*delta_ex) / std::pow(delta_max,2);
  double x1 = 0, x2 = delta_max, x3;
  int count=0;
  double tol = 1.0E-15;
  if (d <= delta_max) {
    while (count < 3000) {
      x3 = (x1+x2)*0.5;
      double a1 = VolumetricHead(x1,a,b,d);
      double a3 = VolumetricHead(x3,a,b,d);

      if (a1*a3 < 0) x2 = x3;
      else x1 = x3;

      if (std::fabs(VolumetricHead(x3,a,b,d)) < tol) break;
      else count++;
    }
  } else {
    return d + delta_ex;
  }
  return x3;
}


double
WeakMPCSemiCoupledDeform::VolumetricHead(double x, double a, double b, double d)
{
  double r = a * std::pow(x,3) + b * std::pow(x,2) - d;
  return r;
}


} // namespace Amanzi


