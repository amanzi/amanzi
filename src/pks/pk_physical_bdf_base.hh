/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Standard base for most PKs, this combines both domains/meshes of
PKPhysicalBase and BDF methods of PKBDFBase.
------------------------------------------------------------------------- */

#ifndef AMANZI_PK_PHYSICAL_BDF_BASE_HH_
#define AMANZI_PK_PHYSICAL_BDF_BASE_HH_

#include "pk_default_base.hh"
#include "pk_bdf_base.hh"
#include "pk_physical_base.hh"

namespace Amanzi {

class PKPhysicalBDFBase : public PKBDFBase, public PKPhysicalBase {

 public:
  PKPhysicalBDFBase(Teuchos::ParameterList& plist,
                    const Teuchos::RCP<TreeVector>& solution) :
      PKDefaultBase(plist, solution),
      PKPhysicalBase(plist,solution),
      PKBDFBase(plist,solution) {}

  virtual void setup(const Teuchos::Ptr<State>& S);

  // initialize.  Note both BDFBase and PhysicalBase have initialize()
  // methods, so we need a unique overrider.
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // Default implementations of BDFFnBase methods.
  // -- Compute a norm on u-du and return the result.
  virtual double enorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du);

  // new virtual set_states() to also get the primary field evaulator.
  virtual void set_states(const Teuchos::RCP<const State>& S,
          const Teuchos::RCP<State>& S_inter,
          const Teuchos::RCP<State>& S_next);

  // -- Experimental approach -- calling this indicates that the time
  //    integration scheme is changing the value of the solution in
  //    state.
  virtual void changed_solution();

 protected:
  // error criteria
  double atol_, rtol_;
  double atol0_, rtol0_;
  bool adapt_tols_to_h_;
  double min_tol_h_;


};


} // namespace

#endif
