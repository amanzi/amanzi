/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Process kernel for energy equation for Richard's flow.
------------------------------------------------------------------------- */

#ifndef PKS_CARBON_SIMPLE_HH_
#define PKS_CARBON_SIMPLE_HH_

#include "PK_Factory.hh"
#include "pk_physical_explicit_default.hh"
#include "PK.hh"

namespace Amanzi {
namespace BGC {

class CarbonSimple : public PK_Physical_Explicit_Default {
 public:
  CarbonSimple(Teuchos::ParameterList& pk_tree,
               const Teuchos::RCP<Teuchos::ParameterList>& glist,
               const Teuchos::RCP<State>& S,
               const Teuchos::RCP<TreeVector>& solution);

  // EnergyBase is a PK
  // -- Setup data
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  // default ok?
  // virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {};

  // -- Calculate any diagnostics prior to doing vis
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S);

  // EnergyBase is a BDFFnBase
  // computes the non-linear functional f = f(t,u,udot)
  virtual void FunctionalTimeDerivative(const double t, const TreeVector& u, TreeVector& f);

  virtual std::string name(){return "carbon simple";};

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
} // namespace Amanzi

  

#endif  
