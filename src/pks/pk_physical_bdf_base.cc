/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Standard base for most PKs, this combines both domains/meshes of
PKPhysicalBase and BDF methods of PKBDFBase.
------------------------------------------------------------------------- */

#include "boost/math/special_functions/fpclassify.hpp"

#include "pk_physical_bdf_base.hh"

namespace Amanzi {


// -----------------------------------------------------------------------------
// Setup
// -----------------------------------------------------------------------------
void PKPhysicalBDFBase::setup(const Teuchos::Ptr<State>& S) {

  // call the meat of the base constructurs via Setup methods
  PKPhysicalBase::setup(S);
  PKBDFBase::setup(S);

  // convergence criteria
  if (conserved_key_.empty()) {
    if (plist_->isParameter("conserved quantity suffix")) {
      Key conserved_default = getKey(domain_, plist_->get<std::string>("conserved quantity suffix"));
      conserved_key_ = plist_->get<std::string>("conserved quantity key", conserved_default);
    } else {
      conserved_key_ = plist_->get<std::string>("conserved quantity key");
    }
  }
  S->RequireField(conserved_key_)->SetMesh(mesh_)
      ->AddComponent("cell",AmanziMesh::CELL,true);
  S->RequireFieldEvaluator(conserved_key_);

  if (cell_vol_key_.empty()) {
    cell_vol_key_ = plist_->get<std::string>("cell volume key",
            getKey(domain_, "cell_volume"));
  }
  S->RequireField(cell_vol_key_)->SetMesh(mesh_)
      ->AddComponent("cell",AmanziMesh::CELL,true);
  S->RequireFieldEvaluator(cell_vol_key_);
  
  atol_ = plist_->get<double>("absolute error tolerance",1.0);
  rtol_ = plist_->get<double>("relative error tolerance",1.0);
  fluxtol_ = plist_->get<double>("flux error tolerance",1.0);
};


// -----------------------------------------------------------------------------
// initialize.  Note both BDFBase and PhysicalBase have initialize()
// methods, so we need a unique overrider.
// -----------------------------------------------------------------------------
void PKPhysicalBDFBase::initialize(const Teuchos::Ptr<State>& S) {
  // Just calls both subclass's initialize.  NOTE - order is important here --
  // PhysicalBase grabs the primary variable and stuffs it into the solution,
  // which must be done prior to BDFBase initializing the timestepper.
  PKPhysicalBase::initialize(S);
  PKBDFBase::initialize(S);
}


// -----------------------------------------------------------------------------
// Default enorm that uses an abs and rel tolerance to monitor convergence.
// -----------------------------------------------------------------------------
double PKPhysicalBDFBase::ErrorNorm(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<const TreeVector> du) {
  S_next_->GetFieldEvaluator(conserved_key_)->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& conserved = *S_next_->GetFieldData(conserved_key_)
      ->ViewComponent("cell",true);
  const Epetra_MultiVector& cv = *S_next_->GetFieldData(cell_vol_key_)
      ->ViewComponent("cell",true);

  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_MEDIUM))
    *vo_->os() << "ENorm (Infnorm) of: " << conserved_key_ << ": " << std::endl;

  Teuchos::RCP<const CompositeVector> dvec = du->Data();
  double h = S_next_->time() - S_inter_->time();

  double enorm_val = 0.0;
  for (CompositeVector::name_iterator comp=dvec->begin();
       comp!=dvec->end(); ++comp) {
    double enorm_comp = 0.0;
    int enorm_loc = -1;
    const Epetra_MultiVector& dvec_v = *dvec->ViewComponent(*comp, false);

    if (*comp == std::string("cell")) {
      // error done relative to extensive, conserved quantity
      int ncells = dvec->size(*comp,false);
      for (unsigned int c=0; c!=ncells; ++c) {
        double enorm_c = std::abs(h * dvec_v[0][c])
            / (atol_*cv[0][c] + rtol_*std::abs(conserved[0][c]));

        if (enorm_c > enorm_comp) {
          enorm_comp = enorm_c;
          enorm_loc = c;
        }
      }
      
    } else if (*comp == std::string("face")) {
      // error in flux -- relative to cell's extensive conserved quantity
      int nfaces = dvec->size(*comp, false);

      for (unsigned int f=0; f!=nfaces; ++f) {
        AmanziMesh::Entity_ID_List cells;
        mesh_->face_get_cells(f, AmanziMesh::OWNED, &cells);
        double cv_min = cells.size() == 1 ? cv[0][cells[0]]
            : std::min(cv[0][cells[0]],cv[0][cells[1]]);
        double conserved_min = cells.size() == 1 ? conserved[0][cells[0]]
            : std::min(conserved[0][cells[0]],conserved[0][cells[1]]);
      
        double enorm_f = fluxtol_ * h * std::abs(dvec_v[0][f])
            / (atol_*cv_min + rtol_*std::abs(conserved_min));
        if (enorm_f > enorm_comp) {
          enorm_comp = enorm_f;
          enorm_loc = f;
        }
      }

    } else {
      double norm;
      dvec_v.Norm2(&norm);
      ASSERT(norm < 1.e-15);
    }

    // Write out Inf norms too.
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      double infnorm(0.);
      dvec_v.NormInf(&infnorm);

      ENorm_t err;
      ENorm_t l_err;
      l_err.value = enorm_comp;
      l_err.gid = dvec_v.Map().GID(enorm_loc);

      int ierr = MPI_Allreduce(&l_err, &err, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
      ASSERT(!ierr);
      *vo_->os() << "  ENorm (" << *comp << ") = " << err.value << "[" << err.gid << "] (" << infnorm << ")" << std::endl;
    }

    enorm_val = std::max(enorm_val, enorm_comp);
  }

  double enorm_val_l = enorm_val;
  int ierr = MPI_Allreduce(&enorm_val_l, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  ASSERT(!ierr);
  return enorm_val;
};


// -----------------------------------------------------------------------------
// Add a boundary marker to owned faces.
// -----------------------------------------------------------------------------
void
PKPhysicalBDFBase::ApplyBoundaryConditions_(const Teuchos::Ptr<CompositeVector>& u) {
  if (u->HasComponent("face")) {
    Epetra_MultiVector& u_f = *u->ViewComponent("face",false);
    unsigned int nfaces = u_f.MyLength();
    for (unsigned int f=0; f!=nfaces; ++f) {
      if (bc_markers_[f] == Operators::OPERATOR_BC_DIRICHLET) {
        u_f[0][f] = bc_values_[f];
      }
    }
  } else if (u->HasComponent("boundary_face")) {
    const Epetra_Map& vandalay_map = mesh_->exterior_face_map(false);
    const Epetra_Map& face_map = mesh_->face_map(false);

    Epetra_MultiVector& u_bf = *u->ViewComponent("boundary_face",false);
    unsigned int nfaces = u_bf.MyLength();
    for (unsigned int bf=0; bf!=nfaces; ++bf) {
      AmanziMesh::Entity_ID f = face_map.LID(vandalay_map.GID(bf));
      if (bc_markers_[f] == Operators::OPERATOR_BC_DIRICHLET) {
        u_bf[0][bf] = bc_values_[f];
      }
    }
  }    
};


double PKPhysicalBDFBase::BoundaryValue(const Teuchos::RCP<const Amanzi::CompositeVector>& solution, int face_id){
  double value=0.;

  if (solution->HasComponent("face")){
    const Epetra_MultiVector& u = *solution -> ViewComponent("face",false);
    value = u[0][face_id];
  }
  else if  (solution->HasComponent("boundary_face")){
    const Epetra_MultiVector& u = *solution -> ViewComponent("boundary_face",false);
    const Epetra_Map& fb_map = mesh_->exterior_face_map(false);
    const Epetra_Map& f_map = mesh_->face_map(false);

    int face_gid = f_map.GID(face_id);
    int face_lbid = fb_map.LID(face_gid);

    value =  u[0][face_lbid];
  }
  else{
    Errors::Message msg("No component is defined for boundary faces\n");
    Exceptions::amanzi_throw(msg);
  }

  return value;

}


  
// -----------------------------------------------------------------------------
// Experimental approach -- calling this indicates that the time
// integration scheme is changing the value of the solution in
// state.
// -----------------------------------------------------------------------------
void PKPhysicalBDFBase::ChangedSolution() {
  solution_evaluator_->SetFieldAsChanged(S_next_.ptr());
};

} // namespace
