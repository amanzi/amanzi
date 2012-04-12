/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Start of a process kernel for the energy equation to be used in thermal
permafrost.  This starts with the simplification that T > T_freezing, limiting
us to the air-water system.
------------------------------------------------------------------------- */

#ifndef PKS_ENERGY_AIRWATERROCK_HH_
#define PKS_ENERGY_AIRWATERROCK_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "state.hh"
#include "tree_vector.hh"
#include "pk_factory.hh"
#include "PK.hh"
#include "bdf_time_integrator.hh"
#include "bdf_fn_base.hh"

#include "twophase_thermal_conductivity.hh"
#include "internal_energy_linear.hh"
#include "internal_energy_gas.hh"
#include "advection.hh"
#include "matrix_mfd.hh"

namespace Amanzi {
namespace Energy {

class TwoPhase : public PK, public BDFFnBase {

public:

  TwoPhase(Teuchos::ParameterList& plist, const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& solution);

  // TwoPhase is a PK
  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::RCP<State>& S);

  // -- transfer operators -- ONLY COPIES POINTERS
  virtual void state_to_solution(const Teuchos::RCP<State>& S,
                                 const Teuchos::RCP<TreeVector>& solution);
  virtual void solution_to_state(const Teuchos::RCP<TreeVector>& solution,
                                 const Teuchos::RCP<State>& S);

  // -- Choose a time step compatible with physics.
  virtual double get_dt();

  // -- Advance from state S0 to state S1 at time S0.time + dt.
  virtual bool advance(double dt);

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S) {}

  // -- Calculate any diagnostics prior to doing vis
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {}


  // TwoPhase is a BDFFnBase
  // computes the non-linear functional f = f(t,u,udot)
  virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f);

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // computes a norm on u-du and returns the result
  virtual double enorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

private:
  // methods for doing science
  void InternalEnergyGas_(const CompositeVector& temp,
                         const CompositeVector& mol_frac_gas,
                         const Teuchos::RCP<CompositeVector>& int_energy_gas);
  void InternalEnergyLiquid_(const CompositeVector& temp,
                            const Teuchos::RCP<CompositeVector>& int_energy_liquid);
  void InternalEnergyRock_(const CompositeVector& temp,
                            const Teuchos::RCP<CompositeVector>& int_energy_rock);

  void EnthalpyLiquid_(const CompositeVector& int_energy_liquid,
          const CompositeVector& pres, const CompositeVector& mol_dens_liq,
          const Teuchos::RCP<CompositeVector>& enthalpy_liq);

  void ThermalConductivity_(const CompositeVector& porosity,
                           const CompositeVector& saturation_liquid,
                           const Teuchos::RCP<CompositeVector>& thermal_conductivity);

  // helper methods for calling the above methods
  void UpdateSecondaryVariables_(const Teuchos::RCP<State>& S);
  void UpdateThermalConductivity_(const Teuchos::RCP<State>& S);
  void UpdateEnthalpyLiquid_(const Teuchos::RCP<State>& S);

  void AddAccumulation_(const Teuchos::RCP<CompositeVector> f);
  void AddAdvection_(const Teuchos::RCP<State> S,
                     const Teuchos::RCP<CompositeVector> f, bool negate);
  void ApplyDiffusion_(const Teuchos::RCP<State> S, const Teuchos::RCP<CompositeVector> f);

  // methods for applying/using the discretization/operators
  void UpdateBoundaryConditions_();
  void ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& temperature);

  // misc setup information
  Teuchos::ParameterList energy_plist_;
  double dt_;

  // boundary conditions
  Teuchos::RCP<BoundaryFunction> bc_temperature_;
  Teuchos::RCP<BoundaryFunction> bc_flux_;
  std::vector<Operators::Matrix_bc> bc_markers_;
  std::vector<double> bc_values_;

  // constitutive relations
  Teuchos::RCP<EnergyRelations::TwophaseThermalConductivity> thermal_conductivity_model_;
  Teuchos::RCP<EnergyRelations::InternalEnergyGas> internal_energy_gas_model_;
  Teuchos::RCP<EnergyRelations::InternalEnergyLinear> internal_energy_liquid_model_;
  Teuchos::RCP<EnergyRelations::InternalEnergyLinear> internal_energy_rock_model_;

  // operators
  Teuchos::RCP<Operators::Advection> advection_;
  Teuchos::RCP<Operators::MatrixMFD> matrix_;
  Teuchos::RCP<Operators::MatrixMFD> preconditioner_;
  std::vector<WhetStone::Tensor> Ke_; // thermal conductivity, needed as tensor for MFD

  // time integration
  Teuchos::RCP<BDFTimeIntegrator> time_stepper_;
  double atol_;
  double rtol_;

  // factory registration
  static RegisteredPKFactory<TwoPhase> reg_;
};

} // namespace Energy
} // namespace Amanzi

#endif
