/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Process kernel for energy equation for Richard's flow.
------------------------------------------------------------------------- */

#ifndef PKS_CARBON_SIMPLE_HH_
#define PKS_CARBON_SIMPLE_HH_

#include "pk_factory.hh"
#include "pk_physical_explicit_base.hh"
#include "PK.hh"

namespace Amanzi {
namespace BGC {

class CarbonSimple : public PKPhysicalExplicitBase {

 public:

  CarbonSimple(const Teuchos::RCP<Teuchos::ParameterList>& plist,
               Teuchos::ParameterList& FElist,
               const Teuchos::RCP<TreeVector>& solution);

  // EnergyBase is a PK
  // -- Setup data
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  // default ok?
  // virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S) {}

  // -- Calculate any diagnostics prior to doing vis
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S);


  // EnergyBase is a BDFFnBase
  // computes the non-linear functional f = f(t,u,udot)
  virtual void Functional(const double t, const TreeVector& u, TreeVector& f);

 protected:

  virtual void ApplyDiffusion_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& g);
  virtual void AddSources_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& g);
  virtual void AddDecomposition_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& g);

  
 protected:
  int npools_;
  
  Key cell_vol_key_;

  bool is_diffusion_;
  Key div_diff_flux_key_;

  bool is_source_;
  Key source_key_;

  bool is_decomp_;
  Key decomp_key_;

 private:
  // factory registration
  static RegisteredPKFactory<CarbonSimple> reg_;
  
};


} // namespace BGC
} // namespace ATS

  

#endif  
