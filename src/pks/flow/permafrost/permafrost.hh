/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  A base three-phase, thermal Richard's equation with water, water vapor, and
  ice for permafrost applications.

  License: BSD
  Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#ifndef PK_FLOW_PERMAFROST_HH_
#define PK_FLOW_PERMAFROST_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "composite_vector.hh"
#include "tree_vector.hh"
#include "state.hh"
#include "matrix_mfd.hh"
#include "upwinding.hh"
#include "boundary_function.hh"

#include "PK.hh"
#include "pk_factory.hh"
#include "bdf_time_integrator.hh"

#include "wrm.hh"
#include "eos.hh"
#include "eos_vapor_in_gas.hh"
#include "pc_ice_water.hh"

#include "richards.hh"

namespace Amanzi {
namespace Flow {

class Permafrost : public Richards {

public:
  Permafrost() {};

  Permafrost(Teuchos::ParameterList& flow_plist, const Teuchos::RCP<State>& S,
           const Teuchos::RCP<TreeVector>& solution);

  // PK methods
  // -- Initialize owned (dependent) variables.
  //    Most methods inherited from Richards
  virtual void initialize(const Teuchos::RCP<State>& S);

  // Permafrost is a BDFFnBase, but inherits most of its methods from Richards
  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

protected:
  void test_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

  // physical methods
  // -- accumulation term
  virtual void AddAccumulation_(const Teuchos::RCP<CompositeVector>& g);

  // -- update secondary variables from primary variables T,p
  virtual void UpdateSecondaryVariables_(const Teuchos::RCP<State>& S);

  // -- Each constitutive relation has a pair of methods, one which pulls data
  //    from a State and calls the second with that data.  This separation is
  //    on purpose, though it is not currently used.  The rationale for this
  //    design choice is that we may eventually move to a Phalanx-like tree
  //    structure for updating the state, and having a single interface entry
  //    point (i.e. a method taking just a state) would likely be useful for
  //    this.

  // Ice EOS
  virtual void UpdateDensityIce_(const Teuchos::RCP<State>& S);
  virtual void DensityIce_(const Teuchos::RCP<State>& S,
                           const CompositeVector& temp,
                           const CompositeVector& pres,
                           const Teuchos::RCP<CompositeVector>& rho_ice,
                           const Teuchos::RCP<CompositeVector>& n_ice);


  // WRM
  virtual void UpdateSaturation_(const Teuchos::RCP<State>& S);

  virtual void FrozenSaturation_(const Teuchos::RCP<State>& S,
          const CompositeVector& temp,
          const CompositeVector& rho_ice,
          const CompositeVector& n_ice,
          const Teuchos::RCP<CompositeVector>& sat_star);

  virtual void UnfrozenSaturation_(const Teuchos::RCP<State>& S,
          const CompositeVector& pres,
          const double& p_atm,
          const Teuchos::RCP<CompositeVector>& sat_star);
  virtual void DUnfrozenSaturationDp_(const Teuchos::RCP<State>& S,
          const CompositeVector& pres,
          const double& p_atm,
          const Teuchos::RCP<CompositeVector>& dsat_star);

  virtual void UpdateRelativePermeability_(const Teuchos::RCP<State>& S);
  virtual void RelativePermeability_(const Teuchos::RCP<State>& S,
          const CompositeVector& sat_liq,
          const Teuchos::RCP<CompositeVector>& rel_perm);

protected:
  Teuchos::RCP<FlowRelations::EOS> eos_ice_;

  // ice-liquid capillary pressure model, may be different in ice wedges?
  // for now, just a single one...
  Teuchos::RCP<FlowRelations::PCIceWater> pc_ice_liq_model_;

  // factory registration
  static RegisteredPKFactory<Permafrost> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
