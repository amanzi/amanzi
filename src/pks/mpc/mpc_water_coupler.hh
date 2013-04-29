/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the derived MPC for water coupling between surface and subsurface.

To be used with either MPCSurfaceSubsurfaceDirichletCoupler or
MPCSurfaceSubsurfaceFluxCoupler.

------------------------------------------------------------------------- */

#ifndef PKS_MPC_WATER_COUPLER_HH_
#define PKS_MPC_WATER_COUPLER_HH_


#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"

#include "tree_vector.hh"
#include "composite_vector.hh"
#include "pk_factory.hh"
#include "pk_default_base.hh"
#include "FieldEvaluator.hh"


namespace Amanzi {

template<class BaseCoupler>
class MPCWaterCoupler : public BaseCoupler, virtual public PKDefaultBase {
 public:

  MPCWaterCoupler(Teuchos::ParameterList& plist,
                  const Teuchos::RCP<TreeVector>& soln);


  // Hackery hook for inheriting MPCs.
  virtual void PreconPostprocess_(Teuchos::RCP<const TreeVector> u,
          Teuchos::RCP<TreeVector> Pu);

  // Given updates to surface cells, calculate updates to surface faces.
  virtual void PreconUpdateSurfaceFaces_(Teuchos::RCP<const TreeVector> u,
          Teuchos::RCP<TreeVector> Pu);

  virtual bool modify_predictor(double h, Teuchos::RCP<TreeVector> up);

 protected:

  bool cap_the_spurt_;
  bool damp_the_spurt_;
  double damp_the_spurt_coef_;
  double face_limiter_;
  double damping_coef_;
  double damping_cutoff_;
  bool modify_predictor_heuristic_;

 private:
  static RegisteredPKFactory< MPCWaterCoupler<BaseCoupler> > reg_;

};


// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
template<class BaseCoupler>
MPCWaterCoupler<BaseCoupler>::MPCWaterCoupler(Teuchos::ParameterList& plist,
        const Teuchos::RCP<TreeVector>& soln) :
    PKDefaultBase(plist, soln),
    BaseCoupler(plist,soln) {

  damping_coef_ = plist.get<double>("damping coefficient", -1.);
  if (damping_coef_ > 0.) {
    damping_cutoff_ = plist.get<double>("damping cutoff", 0.1);
  }

  modify_predictor_heuristic_ =
      plist.get<bool>("modify predictor with heuristic", false);
  face_limiter_ = plist.get<double>("global face limiter", -1);
  cap_the_spurt_ = plist.get<bool>("cap the spurt", false);
  damp_the_spurt_ = plist.isParameter("spurt damping");
  if (damp_the_spurt_) {
    damp_the_spurt_coef_ = plist.get<double>("spurt damping");
  }
}


// -----------------------------------------------------------------------------
// Assorted hacks to try to get water to work
// -----------------------------------------------------------------------------
template<class BaseCoupler>
void MPCWaterCoupler<BaseCoupler>::PreconPostprocess_(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {

  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(this->domain_pk_name_)->data();
  const Epetra_MultiVector& domain_p_f = *S_next_->GetFieldData("pressure")
      ->ViewComponent("face",false);

  // Alterations to PC'd value from limiting and damping
  const double& patm = *S_next_->GetScalarData("atmospheric_pressure");
  Epetra_MultiVector& domain_Pu_f = *domain_Pu->ViewComponent("face",false);
  const Epetra_MultiVector& surf_p_c = *S_next_->GetFieldData("surface_pressure")
      ->ViewComponent("cell",false);

  // global face limiter
  if (face_limiter_ > 0.) {
    int nfaces = domain_Pu_f.MyLength();
    for (int f=0; f!=nfaces; ++f) {
      if (std::abs(domain_Pu_f[0][f]) > face_limiter_) {
        std::cout << "  LIMITING: dp_old = " << domain_Pu_f[0][f];
        domain_Pu_f[0][f] = domain_Pu_f[0][f] > 0. ? face_limiter_ : -face_limiter_;
        std::cout << ", dp_new = " << domain_Pu_f[0][f] << std::endl;

      }
    }
  }

  // Cap surface corrections
  //   In the case that we are starting to infiltrate on dry ground, the low
  //   initial surface rel perm means that to turn on this source, we would
  //   require a huge pressure gradient (to fight the small rel perm).  The
  //   changing rel perm is not represented in the preconditioner, so the
  //   preconditioner tries match the flux by applying a huge gradient,
  //   resulting in a huge update to the surface pressure.  This phenomenon,
  //   known as the spurt, is capped.
  if (cap_the_spurt_) {
    int ncells_surf =
        this->surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (int cs=0; cs!=ncells_surf; ++cs) {
      AmanziMesh::Entity_ID f =
          this->surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);

      double p_old = domain_p_f[0][f];
      double p_new = p_old - domain_Pu_f[0][f];
      if ((p_new > patm) && (p_old < patm - 0.002) && (std::abs(domain_Pu_f[0][f]) > 10000.)) {
        domain_Pu_f[0][f] = p_old - (patm - 0.001);
        std::cout << "  CAPPING: p_old = " << p_old << ", p_new = " << p_new << ", p_capped = " << p_old - domain_Pu_f[0][f] << std::endl;
      }
    }
  }


  if (damp_the_spurt_) {
    int ncells_surf =
        this->surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    int damp = 0;

    for (int cs=0; cs!=ncells_surf; ++cs) {
      AmanziMesh::Entity_ID f =
          this->surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);
      double p_old = domain_p_f[0][f];
      double p_new = p_old - domain_Pu_f[0][f];
      if ((p_new > patm) && (p_old < patm - 0.002) && (std::abs(domain_Pu_f[0][f]) > 10000.)) {
        damp = 1;
      }
    }

    int damp_my = damp;
    this->surf_mesh_->get_comm()->MaxAll(&damp_my, &damp, 1);
    if (damp) {
      std::cout << "  DAMPING THE SPURT!" << std::endl;
      domain_Pu->Scale(damp_the_spurt_coef_);
    }
  }


  // Damp surface corrections
  if (damping_coef_ > 0.) {
    int ncells_surf =
        this->surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (int cs=0; cs!=ncells_surf; ++cs) {
      AmanziMesh::Entity_ID f =
          this->surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);

      double p_old = domain_p_f[0][f];
      double p_new = p_old - domain_Pu_f[0][f];
      if ((p_new > patm) && (std::abs(domain_Pu_f[0][f]) > damping_cutoff_)) {
        domain_Pu_f[0][f] *= damping_coef_;
        std::cout << "  DAMPING: p_old = " << p_old << ", p_new = " << p_new << ", p_damped = " << p_old - domain_Pu_f[0][f] << std::endl;
      }
    }
  }
}


// -----------------------------------------------------------------------------
// Water on the surface is a mixed variable, with height on faces and pressure
// on cells.  Given updates to surface cells, calculate updates to surface
// faces.
// -----------------------------------------------------------------------------
template<class BaseCoupler>
void MPCWaterCoupler<BaseCoupler>::PreconUpdateSurfaceFaces_(
    Teuchos::RCP<const TreeVector> u,
    Teuchos::RCP<TreeVector> Pu) {

  Teuchos::RCP<CompositeVector> surf_Pu = Pu->SubVector(this->surf_pk_name_)->data();
  Epetra_MultiVector& surf_Pu_c = *surf_Pu->ViewComponent("cell",false);
  Teuchos::RCP<const CompositeVector> surf_u = u->SubVector(this->surf_pk_name_)->data();

  // Calculate delta h on the surface
  Teuchos::RCP<CompositeVector> surf_Ph = Teuchos::rcp(new CompositeVector(*surf_Pu));
  surf_Ph->PutScalar(0.);

  // old ponded depth
  *surf_Ph->ViewComponent("cell",false) = *S_next_->GetFieldData("ponded_depth")->ViewComponent("cell",false);

  // new ponded depth
  S_next_->GetFieldData("surface_pressure",this->surf_pk_name_)
      ->ViewComponent("cell",false)->Update(-1., surf_Pu_c, 1.);
  this->surf_pk_->changed_solution();
  S_next_->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S_next_.ptr(), name_);

  // put delta ponded depth into surf_Ph_cell
  surf_Ph->ViewComponent("cell",false)
      ->Update(-1., *S_next_->GetFieldData("ponded_depth")->ViewComponent("cell",false), 1.);

  // update delta faces
  this->surf_preconditioner_->UpdateConsistentFaceCorrection(*surf_u, surf_Ph.ptr());
  *surf_Pu->ViewComponent("face",false) = *surf_Ph->ViewComponent("face",false);

  if (this->surf_c0_ < surf_u->size("cell",false) && this->surf_c1_ < surf_u->size("cell",false)) {
    //  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true)) {
     AmanziMesh::Entity_ID_List fnums1,fnums0;
     std::vector<int> dirs;
     this->surf_mesh_->cell_get_faces_and_dirs(this->surf_c0_, &fnums0, &dirs);
     this->surf_mesh_->cell_get_faces_and_dirs(this->surf_c1_, &fnums1, &dirs);

     *this->out_ << "  PC*u0 in h cell/face: " << (*surf_Ph)("cell",this->surf_c0_) << ", "
           << (*surf_Pu)("face",fnums0[0]) << std::endl;
     *this->out_ << "  PC*u1 in h cell/face: " << (*surf_Ph)("cell",this->surf_c1_) << ", "
           << (*surf_Pu)("face",fnums1[0]) << std::endl;

  }

  // revert solution so we don't break things
  S_next_->GetFieldData("surface_pressure",this->surf_pk_name_)
      ->ViewComponent("cell",false)->Update(1., surf_Pu_c, 1.);
  this->surf_pk_->changed_solution();
}


// -----------------------------------------------------------------------------
// Stop the spurt in a predictor overshoot.
// -----------------------------------------------------------------------------
template<class BaseCoupler>
bool MPCWaterCoupler<BaseCoupler>::modify_predictor(double h,
        Teuchos::RCP<TreeVector> up) {
  Teuchos::OSTab tab = getOSTab();

  // call the BaseCoupler's modify_predictor(), which calls the sub-PKs modify
  // and ensures the surface and subsurface match.
  bool changed = BaseCoupler::modify_predictor(h, up);

  if (modify_predictor_heuristic_) {
    if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true))
      *out_ << " Modifying predictor with water heuristic" << std::endl;

    const Epetra_MultiVector& surf_u_prev_c =
        *S_->GetFieldData("surface_pressure")->ViewComponent("cell",false);
    const double& patm = *S_next_->GetScalarData("atmospheric_pressure");

    Epetra_MultiVector& domain_u_f =
        *up->SubVector(this->domain_pk_name_)->data()->ViewComponent("face",false);
    Epetra_MultiVector& surf_u_c =
        *up->SubVector(this->surf_pk_name_)->data()->ViewComponent("cell",false);

    int ncells = surf_u_c.MyLength();
    for (int c=0; c!=ncells; ++c) {
      int f = this->surf_mesh_->entity_get_parent(AmanziMesh::CELL, c);

      double dp = surf_u_c[0][c] - surf_u_prev_c[0][c];
      double pnew = surf_u_c[0][c] - patm;
      double pold = surf_u_prev_c[0][c] - patm;

      if (pnew > 0) {
        if (dp > pnew) {
          if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true))
            *out_ << "CHANGING (first over?): p = " << surf_u_c[0][c]
                  << " to " << patm + .001 << std::endl;
          surf_u_c[0][c] = patm + .001;
          domain_u_f[0][f] = surf_u_c[0][c];

        } else if (pold > 0 && dp > pold) {
          if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_HIGH, true))
            *out_ << "CHANGING (second over?): p = " << surf_u_c[0][c]
                  << " to " << patm + 2*pold << std::endl;

          surf_u_c[0][c] = patm + 2*pold;
          domain_u_f[0][f] = surf_u_c[0][c];
        }
      }
    }
    changed = true;
  }

  return changed;
}


} // namespace


#endif
