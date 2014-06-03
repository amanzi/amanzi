/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon, Chonggang Xu

   Simple implementation of CLM's Century model for carbon decomposition and a
   simplified 2-PFT (sedge, moss) vegetation model for creating carbon.
   
   ------------------------------------------------------------------------- */

#ifndef PKS_BGC_SIMPLE_HH_
#define PKS_BGC_SIMPLE_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "VerboseObject.hh"
#include "TreeVector.hh"
#include "Explicit_TI_FnBase.hh"

namespace Amanzi {
namespace BGC {

class BGCSimple : public PKPhysicalBase,
                  public Explicit_TI::fnBase<TreeVector> {

 public:
  BGCSimple(const Teuchos::RCP<Teuchos::ParameterList>& plist,
            Teuchos::ParameterList& FElist,
            const Teuchos::RCP<TreeVector>& solution);

  // is a PK
  // -- Setup data
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- provide a timestep size
  virtual double get_dt() {
    return dt_;
  }

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S) {}

  // -- Update diagnostics for vis.
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {}

  // -- advance via one of a few methods -- should get implemented in PKExplicitBase, not here!
  virtual bool advance(double dt);

  // is a Explicit_TI::fnBase
  // -- computes the  functional f = f(t,u) 
  virtual void Functional(const double t, const TreeVector& u, TreeVector& f);
  
};

} // namespace
} // namespace


#endif
