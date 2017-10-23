/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Process kernel for energy equation for Richard's flow.
------------------------------------------------------------------------- */

#ifndef PKS_ENERGY_SURFACE_ICE_HH_
#define PKS_ENERGY_SURFACE_ICE_HH_

#include "PK_Factory.hh"
#include "energy_base.hh"

namespace Amanzi {

// forward declarations
class Function;
namespace Relations { class EOS; }

namespace Energy {

class IEM; 

class EnergySurfaceIce : public EnergyBase {

public:

  EnergySurfaceIce(Teuchos::ParameterList& FElist,
                   const Teuchos::RCP<Teuchos::ParameterList>& plist,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& solution);

  // -- Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S);

protected:
  // -- setup the evaluators
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);

  // -- get enthalpy as a function of Dirichlet boundary data.  Note that this
  //    will get replaced by a better system when we get maps on the boundary
  //    faces.
  //  virtual void ApplyDirichletBCsToEnthalpy_(const Teuchos::Ptr<State>& S);

  virtual void AddSources_(const Teuchos::Ptr<State>& S,
                           const Teuchos::Ptr<CompositeVector>& f);
  virtual void AddSourcesToPrecon_(const Teuchos::Ptr<State>& S, double h);

 protected:
  // models for evaluating enthalpy manually... remove me once boundary faces get in
  Teuchos::RCP<Amanzi::Relations::EOS> eos_liquid_;
  Teuchos::RCP<Energy::IEM> iem_liquid_;

  // simple heat condution term, q = K_s2a * (Tair - Tsurf)
  // air temperature function of time (not space)
  double K_surface_to_air_;

  // flags
  bool standalone_mode_;
  bool is_energy_source_term_;
  bool is_mass_source_term_;
  bool is_air_conductivity_;
  bool mass_source_only_if_unfrozen_;

  // keys
  Key domain_ss_;
  
private:
  // factory registration
  static RegisteredPKFactory<EnergySurfaceIce> reg_;

};

} // namespace Energy
} // namespace Amanzi

#endif
