/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Process kernel for energy equation for Richard's flow.
------------------------------------------------------------------------- */

#ifndef PKS_ENERGY_SURFACE_ICE_HH_
#define PKS_ENERGY_SURFACE_ICE_HH_

#include "pk_factory.hh"
#include "energy_base.hh"

namespace Amanzi {

// forward declarations
namespace Relations { class EOS; }
namespace Energy { namespace EnergyRelations { class IEM; } }

namespace Energy {

class SurfaceIce : public EnergyBase {

public:
  SurfaceIce(Teuchos::ParameterList& plist, const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~SurfaceIce() {}

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S);

protected:
  // -- setup the evaluators
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);

  // -- get enthalpy as a function of Dirichlet boundary data.  Note that this
  //    will get replaced by a better system when we get maps on the boundary
  //    faces.
  virtual void ApplyDirichletBCsToEnthalpy_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& enth);

  // -- Source terms
  virtual void AddSources_(const Teuchos::Ptr<State>& S,
                           const Teuchos::Ptr<CompositeVector>& f);


 protected:
  // models for evaluating enthalpy
  Teuchos::RCP<Relations::EOS> eos_liquid_;
  Teuchos::RCP<EnergyRelations::IEM> iem_liquid_;

  Key energy_source_from_subsurface_key_;
  Key mass_external_source_key_;


private:
  // factory registration
  static RegisteredPKFactory<SurfaceIce> reg_;

};

} // namespace Energy
} // namespace Amanzi

#endif
