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
  atol_ = plist_.get<double>("absolute error tolerance",1.0);
  rtol_ = plist_.get<double>("relative error tolerance",1.0);
  atol0_ = atol_;
  rtol0_ = rtol_;

  // adapt the tolerances to fit the timestep
  adapt_tols_to_h_ = plist_.get<bool>("adapt tolerances to timestep", "false");
  if (adapt_tols_to_h_) {
    min_tol_h_ = plist_.get<double>("cutoff timestep for adaptive tolerance", 100.0);
  }

  // continuation to steady state enables a gradually shrinking atol/rtol
  continuation_to_ss_ = plist_.get<bool>("continuation to steady state", false);

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
double PKPhysicalBDFBase::enorm(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<const TreeVector> du) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = getOSTab();
  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_MEDIUM, true)) {
    *out_ << "ENorm (Infnorm) of: " << name_ << ": ";
  }

  Teuchos::RCP<const CompositeVector> vec = u->data();
  Teuchos::RCP<const CompositeVector> dvec = du->data();

  double enorm_val = 0.0;
  for (CompositeVector::name_iterator comp=vec->begin();
       comp!=vec->end(); ++comp) {
    double enorm_comp = 0.0;
    for (int id=0; id!=vec->size(*comp,false); ++id) {
      double tmp = abs((*dvec)(*comp,id)) / (atol_+rtol_*abs((*vec)(*comp,id)));
      enorm_comp = std::max<double>(enorm_comp, tmp);
    }

    if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_MEDIUM, true)) {
      double buf(0.);
      MPI_Allreduce(&enorm_comp, &buf, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

      double infnorm_comp;
      dvec->ViewComponent(*comp,false)->NormInf(&infnorm_comp);

      *out_ << *comp << " = " << buf << " (" << infnorm_comp << ")  ";
    }
    enorm_val = std::max<double>(enorm_val, enorm_comp);
  }

  if (out_.get() && includesVerbLevel(verbosity_, Teuchos::VERB_MEDIUM, true)) {
    *out_ << std::endl;
  }

#ifdef HAVE_MPI
  double buf = enorm_val;
  MPI_Allreduce(&buf, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
  return enorm_val;
};


// -----------------------------------------------------------------------------
// Experimental approach -- we must pull out S_next_'s solution_evaluator_ to
// stay current for changed_solution()
// -----------------------------------------------------------------------------
void PKPhysicalBDFBase::set_states(const Teuchos::RCP<const State>& S,
        const Teuchos::RCP<State>& S_inter,
        const Teuchos::RCP<State>& S_next) {
  PKDefaultBase::set_states(S, S_inter, S_next);

  Teuchos::RCP<FieldEvaluator> fm = S_next->GetFieldEvaluator(key_);

#if ENABLE_DBC
  solution_evaluator_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  ASSERT(solution_evaluator_ != Teuchos::null);
#else
  solution_evaluator_ = Teuchos::rcp_static_cast<PrimaryVariableFieldEvaluator>(fm);
#endif

  changed_solution();
};


// -----------------------------------------------------------------------------
// Experimental approach -- calling this indicates that the time
// integration scheme is changing the value of the solution in
// state.
// -----------------------------------------------------------------------------
void PKPhysicalBDFBase::changed_solution() {
  solution_evaluator_->SetFieldAsChanged(S_next_.ptr());
};

// -----------------------------------------------------------------------------
// Apply the preconditioner (default application).
// -----------------------------------------------------------------------------
void PKPhysicalBDFBase::precon(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  preconditioner_->ApplyInverse(*u->data(), Pu->data().ptr());
}

} // namespace
