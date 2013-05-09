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
#include "matrix_mfd_surf.hh"

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
  bool damp_and_cap_the_spurt_;
  bool cap_the_cell_spurt_;
  bool make_cells_consistent_;
  double face_limiter_;
  bool modify_predictor_heuristic_;

  Teuchos::RCP<Epetra_MultiVector> surf_pres_prev_;

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

  // predictor modifications
  modify_predictor_heuristic_ =
      plist.get<bool>("modify predictor with heuristic", false);

  // preconditioner modifications
  face_limiter_ = plist.get<double>("global face limiter", -1);
  cap_the_cell_spurt_ = plist.get<bool>("cap the spurt using cells", false);
  cap_the_spurt_ = plist.get<bool>("cap the spurt", false);
  damp_the_spurt_ = plist.get<bool>("damp the spurt", false);
  damp_and_cap_the_spurt_ = plist.get<bool>("damp and cap the spurt", false);
  if (damp_and_cap_the_spurt_) {
    damp_the_spurt_ = true;
    cap_the_spurt_ = true;
  }

  if (face_limiter_ > 0 || cap_the_spurt_ || damp_the_spurt_) {
    make_cells_consistent_ = plist.get<bool>("ensure consistent cells after face updates", false);
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

  // As surface face saturates, its preconditioner is discontinuous at h = 0.
  // Preconditioned updates calculated to go over-saturated using a
  // preconditioner calculated at an under-saturated value overshoot, spurting
  // way too much water onto the surface.  These are assorted methods of
  // dealing with this phenomenon, called "the spurt".

  // Necessary info
  const double& patm = *S_next_->GetScalarData("atmospheric_pressure");
  Epetra_MultiVector& domain_Pu_f = *domain_Pu->ViewComponent("face",false);
  const Epetra_MultiVector& surf_p_c = *S_next_->GetFieldData("surface_pressure")
      ->ViewComponent("cell",false);

  if (surf_pres_prev_ == Teuchos::null)
    surf_pres_prev_ = Teuchos::rcp(new Epetra_MultiVector(surf_p_c));
  *surf_pres_prev_ = *S_next_->GetFieldData("surface_pressure")
      ->ViewComponent("cell",false);


  // Approach 1: global face limiter on the correction size
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

  double damp = 1.;

  // Approach 3: damping of the spurt -- limit the max oversaturated pressure
  //  using a global damping term.
  if (damp_the_spurt_) {
    int ncells_surf =
        this->surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

    for (int cs=0; cs!=ncells_surf; ++cs) {
      AmanziMesh::Entity_ID f =
          this->surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);
      double p_old = domain_p_f[0][f];
      double p_new = p_old - domain_Pu_f[0][f];
      if ((p_new > patm + 100.) && (p_old < patm)) {
        double my_damp = ((patm + 100) - p_old) / (p_new - p_old);
        damp = std::min(damp, my_damp);
      }
    }

    double proc_damp = damp;
    this->surf_mesh_->get_comm()->MinAll(&proc_damp, &damp, 1);
    if (damp < 1.0) {
      std::cout << "  DAMPING THE SPURT!, coef = " << damp << std::endl;
      domain_Pu->Scale(damp);
    }
  }


  // Approach 2: capping of the spurt -- limit the max oversaturated pressure
  //  if coming from undersaturated.
  if (cap_the_spurt_) {
    int ncells_surf =
        this->surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (int cs=0; cs!=ncells_surf; ++cs) {
      AmanziMesh::Entity_ID f =
          this->surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);

      double p_old = domain_p_f[0][f];
      double p_new = p_old - domain_Pu_f[0][f] / damp;
      if ((p_new > patm + 100.) && (p_old < patm)) {
        domain_Pu_f[0][f] = p_old - (patm + 100.);
        std::cout << "  CAPPING: p_old = " << p_old << ", p_new = " << p_new << ", p_capped = " << p_old - domain_Pu_f[0][f] << std::endl;
      }
    }
  }


  // All 3 approaches modify faces.  Cells can be re-back-substituted with
  // these face corrections to get improved cell corrections.
  if (make_cells_consistent_) {
    this->mfd_preconditioner_->UpdateConsistentCellCorrection(
        *u->SubVector(this->domain_pk_name_)->data(),
        Pu->SubVector(this->domain_pk_name_)->data().ptr());
  }

  // Approach 4 works on cells
  if (cap_the_cell_spurt_) {
    int ncapped_l = 0;

    Epetra_MultiVector& domain_Pu_c = *domain_Pu->ViewComponent("cell",false);
    int ncells_surf =
        this->surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    AmanziMesh::Entity_ID_List cells;
    for (int cs=0; cs!=ncells_surf; ++cs) {
      AmanziMesh::Entity_ID f =
          this->surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);
      this->domain_mesh_->face_get_cells(f, AmanziMesh::OWNED, &cells);
      ASSERT(cells.size() == 1);

      double p_old = domain_p_f[0][f];
      double p_new = p_old - domain_Pu_f[0][f];
      if ((p_new > patm + 100.) && (p_old < patm)) {
        ncapped_l++;
        double my_damp = ((patm + 100) - p_old) / (p_new - p_old);
        domain_Pu_c[0][cells[0]] *= my_damp;
        std::cout << "  CAPPING: p_old = " << p_old << ", p_new = " << p_new << ", dp_cell = " << domain_Pu_c[0][cells[0]] << std::endl;
      }
    }

    int ncapped = ncapped_l;
    this->domain_mesh_->get_comm()->SumAll(&ncapped_l, &ncapped, 1);
    if (ncapped > 0) {
      Teuchos::RCP<const CompositeVector> domain_u = u->SubVector(this->domain_pk_name_)->data();
      this->mfd_preconditioner_->UpdateConsistentFaceCorrection(*domain_u, domain_Pu.ptr());
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

  if (this->out_.get() && includesVerbLevel(this->verbosity_, Teuchos::VERB_HIGH, true)) {

    for (std::vector<int>::const_iterator c0=this->surf_dc_.begin();
         c0!=this->surf_dc_.end(); ++c0) {
      if (*c0 < surf_u->size("cell",false)) {
        AmanziMesh::Entity_ID_List fnums0;
        std::vector<int> dirs;
        this->surf_mesh_->cell_get_faces_and_dirs(*c0, &fnums0, &dirs);

        *this->out_ << "  PC*u(" << *c0 << ") in h cell/face: "
                    << (*surf_Ph)("cell",*c0) << ", "
                    << (*surf_Pu)("face",fnums0[0]) << std::endl;
      }
    }
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
