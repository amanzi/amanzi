// Delegate for heuristic corrections based upon coupled surface/subsurface water.

#include "mpc_delegate_water.hh"

namespace Amanzi {

  MPCDelegateWater::MPCDelegateWater(const Teuchos::RCP<Teuchos::ParameterList>& plist, std::string domain) :
    plist_(plist),
    i_domain_(-1),
    i_surf_(-1),
    i_Tdomain_(-1),
    i_Tsurf_(-1)
{
  // predictor control
  modify_predictor_heuristic_ =
      plist_->get<bool>("modify predictor with heuristic", false);
  modify_predictor_spurt_damping_ =
      plist_->get<bool>("modify predictor damp and cap the water spurt", false);
  modify_predictor_tempfromsource_ =
      plist_->get<bool>("modify predictor surface temperature from source", false);

  // precon control
  face_limiter_ = plist_->get<double>("global water face limiter", -1.0);
  cap_the_spurt_ = plist_->get<bool>("cap the water spurt", false);
  damp_the_spurt_ = plist_->get<bool>("damp the water spurt", false);
  bool damp_and_cap_the_spurt = plist_->get<bool>("damp and cap the water spurt", false);
  if (damp_and_cap_the_spurt) {
    damp_the_spurt_ = true;
    cap_the_spurt_ = true;
  }

  cap_the_sat_spurt_ = plist_->get<bool>("cap the saturated spurt", false);
  damp_the_sat_spurt_ = plist_->get<bool>("damp the saturated spurt", false);
  bool damp_and_cap_the_sat_spurt = plist_->get<bool>("damp and cap the saturated spurt", false);
  if (damp_and_cap_the_sat_spurt) {
    damp_the_sat_spurt_ = true;
    cap_the_sat_spurt_ = true;
  }
  
  // set the size of the caps
  if (cap_the_spurt_ || damp_the_spurt_ ||
      cap_the_sat_spurt_ || damp_the_sat_spurt_ ||
      modify_predictor_heuristic_ || modify_predictor_spurt_damping_) {
    cap_size_ = plist_->get<double>("cap over atmospheric", 100.0);
  }

  // create the VO
  vo_ = Teuchos::rcp(new VerboseObject(plist->name(), *plist_));

  domain_ss = domain;
}

// Approach 1: global face limiter on the correction size
int
MPCDelegateWater::ModifyCorrection_WaterFaceLimiter(double h, Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  int n_modified = 0;
  if (face_limiter_ > 0.) {
    Epetra_MultiVector& domain_Pu_f = *Pu->SubVector(i_domain_)->Data()
        ->ViewComponent("face",false);

    for (int f=0; f!=domain_Pu_f.MyLength(); ++f) {
      if (std::abs(domain_Pu_f[0][f]) > face_limiter_) {
        if (vo_->os_OK(Teuchos::VERB_HIGH))
          *vo_->os() << "  LIMITING: dp_old = " << domain_Pu_f[0][f];
        domain_Pu_f[0][f] = domain_Pu_f[0][f] > 0. ? face_limiter_ : -face_limiter_;
        if (vo_->os_OK(Teuchos::VERB_HIGH))
          *vo_->os() << ", dp_new = " << domain_Pu_f[0][f] << std::endl;
        n_modified++;
      }
    }
  }

  return n_modified;
}

// Approach 2: damping of the spurt -- limit the max oversaturated pressure
//  using a global damping term.
// Approach 3: capping of the spurt -- limit the max oversaturated pressure
//  if coming from undersaturated.
double
MPCDelegateWater::ModifyCorrection_WaterSpurtDamp(double h, Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  const double& patm = *S_next_->GetScalarData("atmospheric_pressure");

  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh =
      u->SubVector(i_surf_)->Data()->Mesh();
  int ncells_surf = surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  const Epetra_MultiVector& domain_p_f = *u->SubVector(i_domain_)->Data()
      ->ViewComponent("face",false);
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(i_domain_)->Data();
  Epetra_MultiVector& domain_Pu_f = *domain_Pu->ViewComponent("face",false);

  // Approach 2
  double damp = 1.;
  if (damp_the_spurt_) {
    for (int cs=0; cs!=ncells_surf; ++cs) {
      AmanziMesh::Entity_ID f =
          surf_mesh->entity_get_parent(AmanziMesh::CELL, cs);
      double p_old = domain_p_f[0][f];
      double p_new = p_old - domain_Pu_f[0][f];
      if ((p_new > patm + cap_size_) && (p_old < patm)) {
        double my_damp = ((patm + cap_size_) - p_old) / (p_new - p_old);
        damp = std::min(damp, my_damp);
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          std::cout << "   DAMPING THE SPURT (sc=" << surf_mesh->cell_map(false).GID(cs) << "): p_old = " << p_old << ", p_new = " << p_new << ", coef = " << my_damp << std::endl;
      }
    }

    double proc_damp = damp;
    domain_Pu_f.Comm().MinAll(&proc_damp, &damp, 1);
    if (damp < 1.0) {
      if (vo_->os_OK(Teuchos::VERB_HIGH))
        *vo_->os() << "  DAMPING THE SPURT!, coef = " << damp << std::endl;
      domain_Pu->Scale(damp);
    }
  }
  return damp;
}


int
MPCDelegateWater::ModifyCorrection_WaterSpurtCap(double h, Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu, double damp) {
  const double& patm = *S_next_->GetScalarData("atmospheric_pressure");

  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh =
      u->SubVector(i_surf_)->Data()->Mesh();
  int ncells_surf = surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  const Epetra_MultiVector& domain_p_f = *u->SubVector(i_domain_)->Data()
      ->ViewComponent("face",false);
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(i_domain_)->Data();
  Epetra_MultiVector& domain_Pu_f = *domain_Pu->ViewComponent("face",false);

  // Approach 3
  int n_modified = 0;
  if (cap_the_spurt_) {
    for (int cs=0; cs!=ncells_surf; ++cs) {
      AmanziMesh::Entity_ID f = surf_mesh->entity_get_parent(AmanziMesh::CELL, cs);

      double p_old = domain_p_f[0][f];
      double p_new = p_old - domain_Pu_f[0][f] / damp;
      if ((p_new > patm + cap_size_) && (p_old < patm)) {
        domain_Pu_f[0][f] = p_old - (patm + cap_size_);
        n_modified++;
        if (vo_->os_OK(Teuchos::VERB_HIGH))
          std::cout << "  CAPPING THE SPURT (sc=" << surf_mesh->cell_map(false).GID(cs) << ",f=" 
		    << u->SubVector(i_domain_)->Data()->Mesh()->face_map(false).GID(f) << "): p_old = " << p_old
                    << ", p_new = " << p_new << ", p_capped = " << p_old - domain_Pu_f[0][f] << std::endl;
      } else if ((p_new < patm) && (p_old > patm)) {
        // strange attempt to kick NKA when it goes back under?
        n_modified++;
        if (vo_->os_OK(Teuchos::VERB_HIGH))
          std::cout << "  INVERSE SPURT (sc=" << surf_mesh->cell_map(false).GID(cs) << "): p_old = " << p_old
                    << ", p_new = " << p_new << std::endl;
      }
    }
  }

  return n_modified;
}


// Approach 2: damping of the spurt -- limit the max oversaturated pressure
//  using a global damping term.
// Approach 3: capping of the spurt -- limit the max oversaturated pressure
//  if coming from undersaturated.
double
MPCDelegateWater::ModifyCorrection_SaturatedSpurtDamp(double h, Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  const double& patm = *S_next_->GetScalarData("atmospheric_pressure");

  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh =
      u->SubVector(i_domain_)->Data()->Mesh();
  int ncells_domain = domain_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  const Epetra_MultiVector& domain_p_c = *u->SubVector(i_domain_)->Data()
      ->ViewComponent("cell",false);
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(i_domain_)->Data();
  Epetra_MultiVector& domain_Pu_c = *domain_Pu->ViewComponent("cell",false);

  // Approach 2
  double damp = 1.;
  if (damp_the_sat_spurt_) {
    for (int c=0; c!=ncells_domain; ++c) {
      double p_old = domain_p_c[0][c];
      double p_new = p_old - domain_Pu_c[0][c];
      if ((p_new > patm + cap_size_) && (p_old < patm)) {
        double my_damp = ((patm + cap_size_) - p_old) / (p_new - p_old);
        damp = std::min(damp, my_damp);
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          std::cout << "   DAMPING THE SATURATED SPURT (c=" << domain_mesh->cell_map(false).GID(c) << "): p_old = " << p_old << ", p_new = " << p_new << ", coef = " << my_damp << std::endl;
      }
    }

    double proc_damp = damp;
    domain_Pu_c.Comm().MinAll(&proc_damp, &damp, 1);
    if (damp < 1.0) {
      if (vo_->os_OK(Teuchos::VERB_HIGH))
        *vo_->os() << "  DAMPING THE SATURATED SPURT!, coef = " << damp << std::endl;
      domain_Pu->Scale(damp);
    }
  }
  return damp;
}


int
MPCDelegateWater::ModifyCorrection_SaturatedSpurtCap(double h, Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu, double damp) {
  const double& patm = *S_next_->GetScalarData("atmospheric_pressure");

  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh =
      u->SubVector(i_domain_)->Data()->Mesh();
  int ncells_domain = domain_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  const Epetra_MultiVector& domain_p_c = *u->SubVector(i_domain_)->Data()
      ->ViewComponent("cell",false);
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(i_domain_)->Data();
  Epetra_MultiVector& domain_Pu_c = *domain_Pu->ViewComponent("cell",false);

  // Approach 3
  int n_modified = 0;
  if (cap_the_sat_spurt_) {
    for (int c=0; c!=ncells_domain; ++c) {
      double p_old = domain_p_c[0][c];
      double p_new = p_old - domain_Pu_c[0][c] / damp;
      if ((p_new > patm + cap_size_) && (p_old < patm)) {
        domain_Pu_c[0][c] = p_old - (patm + cap_size_);
        n_modified++;
        if (vo_->os_OK(Teuchos::VERB_HIGH))
          std::cout << "  CAPPING THE SATURATED SPURT (c=" << domain_mesh->cell_map(false).GID(c)
                    << "): p_old = " << p_old
                    << ", p_new = " << p_new << ", p_capped = " << p_old - domain_Pu_c[0][c] << std::endl;
      }
    }
  }

  return n_modified;
}

// modify predictor via heuristic stops spurting in the surface flow
bool
MPCDelegateWater::ModifyPredictor_Heuristic(double h, const Teuchos::RCP<TreeVector>& u) {
  bool modified = false;
  if (modify_predictor_heuristic_) {
    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os() << "  MPCWaterCoupler: Modifying predictor with water heuristic" << std::endl;

    Epetra_MultiVector& domain_u_f =
        *u->SubVector(i_domain_)->Data()->ViewComponent("face",false);
    Epetra_MultiVector& surf_u_c =
        *u->SubVector(i_surf_)->Data()->ViewComponent("cell",false);

    Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh =
        u->SubVector(i_surf_)->Data()->Mesh();

    const Epetra_MultiVector& surf_u_prev_c =
        *S_->GetFieldData("surface_pressure")->ViewComponent("cell",false);
    const double& patm = *S_next_->GetScalarData("atmospheric_pressure");
    int ncells = surf_u_c.MyLength();
    for (int c=0; c!=ncells; ++c) {
      int f = surf_mesh->entity_get_parent(AmanziMesh::CELL, c);

      double dp = surf_u_c[0][c] - surf_u_prev_c[0][c];
      double pnew = surf_u_c[0][c] - patm;
      double pold = surf_u_prev_c[0][c] - patm;

      if (pnew > 0) {
        if (dp > pnew) {
          if (vo_->os_OK(Teuchos::VERB_HIGH))
            *vo_->os() << "CHANGING (first over?): p = " << surf_u_c[0][c]
                       << " to " << patm + cap_size_ << std::endl;
          surf_u_c[0][c] = patm + cap_size_;
          domain_u_f[0][f] = surf_u_c[0][c];

        } else if (pold > 0 && dp > pold) {
          if (vo_->os_OK(Teuchos::VERB_HIGH))
            *vo_->os() << "CHANGING (second over?): p = " << surf_u_c[0][c]
                       << " to " << patm + 2*pold << std::endl;
          surf_u_c[0][c] = patm + 2*pold;
          domain_u_f[0][f] = surf_u_c[0][c];
        }
      }
    }
    modified = true;
  }

  return modified;
}

// Approach 2: damping of the spurt -- limit the max oversaturated pressure
//  using a global damping term.
// Approach 3: capping of the spurt -- limit the max oversaturated pressure
//  if coming from undersaturated.
bool
MPCDelegateWater::ModifyPredictor_WaterSpurtDamp(double h,
        const Teuchos::RCP<TreeVector>& u) {

  // Approach 2
  if (modify_predictor_spurt_damping_) {
    const double& patm = *S_next_->GetScalarData("atmospheric_pressure");

    Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh =
        u->SubVector(i_surf_)->Data()->Mesh();
    int ncells_surf = surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

    Epetra_MultiVector& surf_pnew_c = *u->SubVector(i_surf_)->Data()
        ->ViewComponent("cell",false);
    Epetra_MultiVector& domain_pnew_f = *u->SubVector(i_domain_)->Data()
        ->ViewComponent("face",false);
    Key key_ss;
    if (domain_ss.substr(0,6) == "column")
      key_ss = getKey(domain_ss,"pressure");
    else
      key_ss = "pressure";

    const Epetra_MultiVector& domain_pold_f =
        *S_->GetFieldData(key_ss)->ViewComponent("face",false);

    int rank = surf_mesh->get_comm()->MyPID();
    double damp = 1.;
    for (unsigned int cs=0; cs!=ncells_surf; ++cs) {
      AmanziMesh::Entity_ID f =
          surf_mesh->entity_get_parent(AmanziMesh::CELL, cs);
      double p_old = domain_pold_f[0][f];
      double p_new = domain_pnew_f[0][f];
      if ((p_new > patm + cap_size_) && (p_old < patm)) {
        // first over
        double my_damp = ((patm + cap_size_) - p_old) / (p_new - p_old);
        damp = std::min(damp, my_damp);
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          std::cout << "   DAMPING THE SPURT (1st over) (sc=" << surf_mesh->cell_map(false).GID(cs) << "): p_old = " << p_old << ", p_new = " << p_new << ", coef = " << my_damp << std::endl;
      } else if ((p_old > patm) && (p_new - p_old > p_old - patm)) {
        // second over
        double my_damp = ((patm + 2*(p_old - patm)) - p_old) / (p_new - p_old);
        damp = std::min(damp, my_damp);
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          std::cout << "   DAMPING THE SPURT (2nd over) (sc=" << surf_mesh->cell_map(false).GID(cs) << "): p_old = " << p_old << ", p_new = " << p_new << ", coef = " << my_damp << std::endl;
      }
    }

    double proc_damp = damp;
    domain_pnew_f.Comm().MinAll(&proc_damp, &damp, 1);
    if (damp < 1.0) {
      if (vo_->os_OK(Teuchos::VERB_HIGH))
        *vo_->os() << "  DAMPING THE SPURT!, coef = " << damp << std::endl;

      // apply the damping
      Teuchos::RCP<const CompositeVector> domain_pold = S_->GetFieldData(key_ss);
      Teuchos::RCP<CompositeVector> domain_pnew = u->SubVector(i_domain_)->Data();
      db_->WriteVector("p_old", domain_pold.ptr());
      db_->WriteVector("p_new", domain_pnew.ptr());
      domain_pnew->Update(1. - damp, *domain_pold, damp);
      db_->WriteVector("p_damped", domain_pnew.ptr());

      // undamp and cap the surface
      for (unsigned int cs=0; cs!=ncells_surf; ++cs) {
        AmanziMesh::Entity_ID f =
            surf_mesh->entity_get_parent(AmanziMesh::CELL, cs);
        double p_old = domain_pold_f[0][f];
        double p_new = (domain_pnew_f[0][f] - p_old) / damp + p_old;
        if ((p_new > patm + cap_size_) && (p_old < patm)) {
          // first over
          domain_pnew_f[0][f] = patm + cap_size_;
          surf_pnew_c[0][cs] = domain_pnew_f[0][f];
          if (vo_->os_OK(Teuchos::VERB_HIGH))
            std::cout << "  CAPPING THE SPURT (1st over) (sc=" << surf_mesh->cell_map(false).GID(cs) << "): p_old = " << p_old << ", p_new = " << p_new << ", p_capped = " << domain_pnew_f[0][f] << std::endl;
        } else if ((p_old > patm) && (p_new - p_old > p_old - patm)) {
          // second over
          domain_pnew_f[0][f] = patm + 2*(p_old - patm);
          surf_pnew_c[0][cs] = domain_pnew_f[0][f];
          if (vo_->os_OK(Teuchos::VERB_HIGH))
            std::cout << "  CAPPING THE SPURT (2nd over) (sc=" << surf_mesh->cell_map(false).GID(cs) << "): p_old = " << p_old << ", p_new = " << p_new << ", p_capped = " << domain_pnew_f[0][f] << std::endl;
        } else {
          surf_pnew_c[0][cs] = domain_pnew_f[0][f];
        }
      }
    }

    return damp < 1.0;
  }

  return false;
}


// modify predictor via heuristic stops spurting in the surface flow
bool
MPCDelegateWater::ModifyPredictor_TempFromSource(double h, const Teuchos::RCP<TreeVector>& u) {
  bool modified = false;
  if (modify_predictor_tempfromsource_) {
    modified = true;
    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os() << "  MPCWaterCoupler: Modifying predictor, taking surface temperature from source." << std::endl;

    Epetra_MultiVector& surf_Tnew_c = *u->SubVector(i_Tsurf_)->Data()
        ->ViewComponent("cell",false);
    Epetra_MultiVector& surf_pnew_c = *u->SubVector(i_surf_)->Data()
        ->ViewComponent("cell",false);
    Epetra_MultiVector& domain_pnew_f = *u->SubVector(i_domain_)->Data()
        ->ViewComponent("face",false);
    const Epetra_MultiVector& Told = *S_->GetFieldData("surface_temperature")
        ->ViewComponent("cell",false);
    const Epetra_MultiVector& Tsource = *S_next_->GetFieldData("surface_mass_source_temperature")
        ->ViewComponent("cell",false);
    const Epetra_MultiVector& hold = *S_->GetFieldData("ponded_depth")
        ->ViewComponent("cell",false);
    const Epetra_MultiVector& dhsource = *S_->GetFieldData("surface_mass_source")
        ->ViewComponent("cell",false);

    Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh =
        u->SubVector(i_surf_)->Data()->Mesh();

    for (unsigned int c=0; c!=surf_Tnew_c.MyLength(); ++c) {
      if (surf_Tnew_c[0][c] < 271.15) {
        // frozen, modify predictor to ensure surface is ready to accept ice
        if (surf_pnew_c[0][c] < 101325.) {
          surf_pnew_c[0][c] = 101325.1;
          AmanziMesh::Entity_ID f =
              surf_mesh->entity_get_parent(AmanziMesh::CELL, c);
          domain_pnew_f[0][f] = surf_pnew_c[0][c];
        }
      }
    }
  }
  return modified;
}

} // namespace
