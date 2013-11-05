// Delegate for heuristic corrections based upon coupled surface/subsurface water.

#include "mpc_delegate_water.hh"

namespace Amanzi {

MPCDelegateWater::MPCDelegateWater(const Teuchos::RCP<Teuchos::ParameterList>& plist,
        int i_domain, int i_surf) :
    plist_(plist),
    i_domain_(i_domain),
    i_surf_(i_surf)
{

  face_limiter_ = plist_->get<double>("global water face limiter", -1.0);
  cap_the_spurt_ = plist_->get<bool>("cap the water spurt", false);
  damp_the_spurt_ = plist_->get<bool>("damp the water spurt", false);
  bool damp_and_cap_the_spurt = plist_->get<bool>("damp and cap the water spurt", false);
  if (damp_and_cap_the_spurt) {
    damp_the_spurt_ = true;
    cap_the_spurt_ = true;
  }

  if (cap_the_spurt_ || damp_the_spurt_) {
    cap_size_ = plist_->get<double>("cap over atmospheric", 100.0);
  }

  vo_ = Teuchos::rcp(new VerboseObject(plist->name(), *plist_));
}

bool
MPCDelegateWater::ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du) {
  bool modified = false;

  int n_modified = 0;
  if (face_limiter_ > 0.) {
    n_modified += ModifyCorrection_WaterFaceLimiter_(h, res, u, du);
  }

  if (cap_the_spurt_ || damp_the_spurt_) {
    n_modified += ModifyCorrection_WaterSpurt_(h, res, u, du);
  }

  if (cap_the_spurt_ || damp_the_spurt_ || face_limiter_ > 0.) {
    int n_modified_l = n_modified;
    u->SubVector(i_domain_)->Data()->Comm().SumAll(&n_modified_l, &n_modified, 1);
    if (n_modified > 0) modified = true;
  }

  return modified;
}

// Approach 1: global face limiter on the correction size
int
MPCDelegateWater::ModifyCorrection_WaterFaceLimiter_(double h, Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  int n_modified = 0;
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

  return n_modified;
}


// Approach 2: damping of the spurt -- limit the max oversaturated pressure
//  using a global damping term.
// Approach 3: capping of the spurt -- limit the max oversaturated pressure
//  if coming from undersaturated.
int
MPCDelegateWater::ModifyCorrection_WaterSpurt_(double h, Teuchos::RCP<const TreeVector> res,
        Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) {
  double patm = 101325.;

  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh =
      u->SubVector(i_surf_)->Data()->Mesh();
  int ncells_surf = surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  const Epetra_MultiVector& domain_p_f = *u->SubVector(i_domain_)->Data()
      ->ViewComponent("face",false);
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(i_domain_)->Data();
  Epetra_MultiVector& domain_Pu_f = *domain_Pu->ViewComponent("face",false);

  // Approach 2
  int n_affected = 0;
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
      }
    }

    double proc_damp = damp;
    domain_Pu_f.Comm().MinAll(&proc_damp, &damp, 1);
    if (damp < 1.0) {
      if (vo_->os_OK(Teuchos::VERB_HIGH))
        *vo_->os() << "  DAMPING THE SPURT!, coef = " << damp << std::endl;
      domain_Pu->Scale(damp);
      n_affected = 1;
    }
  }

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
          *vo_->os() << "  CAPPING THE SPURT: p_old = " << p_old << ", p_new = " << p_new << ", p_capped = " << p_old - domain_Pu_f[0][f] << std::endl;
      }
    }
  }

  return n_affected + n_modified;
}


} // namespace
