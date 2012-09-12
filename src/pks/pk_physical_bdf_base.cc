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
  PKBDFBase::setup(S);
  PKPhysicalBase::setup(S);
};


// -----------------------------------------------------------------------------
// initialize.  Note both BDFBase and PhysicalBase have initialize()
// methods, so we need a unique overrider.
// -----------------------------------------------------------------------------
void PKPhysicalBDFBase::initialize(const Teuchos::Ptr<State>& S) {
  // just calls both subclasses
  PKBDFBase::initialize(S);
  PKPhysicalBase::initialize(S);
}


// -----------------------------------------------------------------------------
// Default enorm that uses an abs and rel tolerance to monitor convergence.
// -----------------------------------------------------------------------------
double PKPhysicalBDFBase::enorm(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<const TreeVector> du) {
  Teuchos::RCP<const CompositeVector> vec = u->data();
  Teuchos::RCP<const CompositeVector> dvec = du->data();

  double enorm_val = 0.0;
  for (CompositeVector::name_iterator comp=vec->begin();
       comp!=vec->end(); ++comp) {
    for (int id=0; id!=vec->size(*comp,false); ++id) {
      if (boost::math::isnan<double>((*dvec)(*comp,id))) {
        std::cout << "Cutting time step due to NaN in correction." << std::endl;
        Errors::Message m("Cut time step");
        Exceptions::amanzi_throw(m);
      }

      double tmp = abs((*dvec)(*comp,id)) / (atol_+rtol_*abs((*vec)(*comp,id)));
      enorm_val = std::max<double>(enorm_val, tmp);
    }
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
};


// -----------------------------------------------------------------------------
// Experimental approach -- calling this indicates that the time
// integration scheme is changing the value of the solution in
// state.
// -----------------------------------------------------------------------------
void PKPhysicalBDFBase::changed_solution() {
  solution_evaluator_->SetFieldAsChanged();
};

} // namespace
