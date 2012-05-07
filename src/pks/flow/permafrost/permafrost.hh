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

#include "PK.hh"
#include "pk_factory.hh"
#include "bdf_time_integrator.hh"

#include "wrm.hh"
#include "eos.hh"
#include "eos_vapor_in_gas.hh"
#include "pc_ice_water.hh"

namespace Amanzi {
namespace Flow {

const int FLOW_RELATIVE_PERM_CENTERED = 1;
const int FLOW_RELATIVE_PERM_UPWIND_GRAVITY = 2;
const int FLOW_RELATIVE_PERM_UPWIND_DARCY_FLUX = 3;
const int FLOW_RELATIVE_PERM_ARITHMETIC_MEAN = 4;

class Permafrost : public PK {

public:
  Permafrost(Teuchos::ParameterList& flow_plist, const Teuchos::RCP<State>& S,
           const Teuchos::RCP<TreeVector>& solution);

  // main methods
  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::RCP<State>& S);

  // -- Choose a time step compatible with physics.
  virtual double get_dt() {
    return dt_;
  }

  // -- transfer operators -- pointer copy
  virtual void state_to_solution(const Teuchos::RCP<State>& S,
          const Teuchos::RCP<TreeVector>& soln);
  virtual void solution_to_state(const Teuchos::RCP<TreeVector>& soln,
          const Teuchos::RCP<State>& S);

  // -- Advance from state S to state S_next at time S0.time + dt.
  virtual bool advance(double dt);

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S);

  // -- Update diagnostics for vis.
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S);

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // computes a norm on u-du and returns the result
  virtual double enorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

private:
  // boundary condition members
  virtual void UpdateBoundaryConditions_();
  virtual void ApplyBoundaryConditions_(const Teuchos::RCP<State>& S,
          const Teuchos::RCP<CompositeVector>& temperature);

  // bdf needs help
  void DeriveFaceValuesFromCellValues_(const Teuchos::RCP<State>& S,
                                       const Teuchos::RCP<CompositeVector>& pres);

  // computational concerns in managing abs, rel perm
  // -- is abs perm changing?
  bool variable_abs_perm() { return variable_abs_perm_; }

  // -- builds tensor K, along with faced-based Krel if needed by the rel-perm method
  void SetAbsolutePermeabilityTensor_(const Teuchos::RCP<State>& S);
  void UpdatePermeabilityData_(const Teuchos::RCP<State>& S);

  // -- rel perm calculation for fluxes
  void CalculateRelativePermeabilityUpwindGravity_(const Teuchos::RCP<State>& S,
          const CompositeVector& rel_perm_cells,
          const Teuchos::RCP<CompositeVector>& rel_perm_faces);
  void CalculateRelativePermeabilityUpwindFlux_(const Teuchos::RCP<State>& S,
          const CompositeVector& flux, const CompositeVector& rel_perm_cells,
          const Teuchos::RCP<CompositeVector>& rel_perm_faces);
  void CalculateRelativePermeabilityArithmeticMean_(const Teuchos::RCP<State>& S,
          const CompositeVector& rel_perm_cells,
          const Teuchos::RCP<CompositeVector>& rel_perm_faces);

  // physical methods
  // -- diffusion term
  void ApplyDiffusion_(const Teuchos::RCP<State>& S,const Teuchos::RCP<CompositeVector>& g);

  // -- accumulation term
  void AddAccumulation_(const Teuchos::RCP<CompositeVector>& g);

  // -- gravity contributions
  void AddGravityFluxes_(const Teuchos::RCP<State>& S,
          const Teuchos::RCP<Operators::MatrixMFD>& matrix);
  void AddGravityFluxesToVector_(const Teuchos::RCP<State>& S,
          const Teuchos::RCP<CompositeVector>& darcy_mass_flux);

  // -- update secondary variables from primary variables T,p
  void UpdateSecondaryVariables_(const Teuchos::RCP<State>& S);

  // -- Each constitutive relation has a pair of methods, one which pulls data
  //    from a State and calls the second with that data.  This separation is
  //    on purpose, though it is not currently used.  The rationale for this
  //    design choice is that we may eventually move to a Phalanx-like tree
  //    structure for updating the state, and having a single interface entry
  //    point (i.e. a method taking just a state) would likely be useful for
  //    this.

  // Liquid EOS
  void UpdateDensityLiquid_(const Teuchos::RCP<State>& S);
  void DensityLiquid_(const Teuchos::RCP<State>& S,
                      const CompositeVector& temp,
                      const CompositeVector& pres,
                      const Teuchos::RCP<CompositeVector>& rho_liq,
                      const Teuchos::RCP<CompositeVector>& n_liq);

  void UpdateViscosityLiquid_(const Teuchos::RCP<State>& S);
  void ViscosityLiquid_(const Teuchos::RCP<State>& S,
                        const CompositeVector& temp,
                        const Teuchos::RCP<CompositeVector>& visc_liq);

  // Ice EOS
  void UpdateDensityIce_(const Teuchos::RCP<State>& S);
  void DensityIce_(const Teuchos::RCP<State>& S,
                   const CompositeVector& temp,
                   const CompositeVector& pres,
                   const Teuchos::RCP<CompositeVector>& rho_ice,
                   const Teuchos::RCP<CompositeVector>& n_ice);

  // Vapor EOS
  void UpdateDensityGas_(const Teuchos::RCP<State>& S);
  void DensityGas_(const Teuchos::RCP<State>& S,
                   const CompositeVector& temp,
                   const CompositeVector& pres,
                   const double& p_atm,
                   const Teuchos::RCP<CompositeVector>& mol_frac_gas,
                   const Teuchos::RCP<CompositeVector>& rho_gas,
                   const Teuchos::RCP<CompositeVector>& n_gas);

  // WRM
  void UpdateSaturation_(const Teuchos::RCP<State>& S);

  void FrozenSaturation_(const Teuchos::RCP<State>& S,
                         const CompositeVector& temp,
                         const CompositeVector& rho_ice,
                         const CompositeVector& n_ice,
                         const Teuchos::RCP<CompositeVector>& sat_star);

  void UnfrozenSaturation_(const Teuchos::RCP<State>& S,
                           const CompositeVector& pres,
                           const double& p_atm,
                           const Teuchos::RCP<CompositeVector>& sat_star);
  void DUnfrozenSaturationDp_(const Teuchos::RCP<State>& S,
          const CompositeVector& pres,
          const double& p_atm,
          const Teuchos::RCP<CompositeVector>& dsat_star);

  void UpdateRelativePermeability_(const Teuchos::RCP<State>& S);
  void RelativePermeability_(const Teuchos::RCP<State>& S,
                             const CompositeVector& sat_liq,
                             const Teuchos::RCP<CompositeVector>& rel_perm);

private:
  // control switches
  int internal_tests_;  // output information
  bool variable_abs_perm_;
  int Krel_method_;

  double dt_;
  double dt0_;

  // input parameter data
  Teuchos::ParameterList flow_plist_;

  // work data space
  std::vector<WhetStone::Tensor> K_;  // tensor of absolute permeability

  // constitutive relations
  Teuchos::RCP<FlowRelations::EOS> eos_liquid_;
  Teuchos::RCP<FlowRelations::EOS> eos_ice_;
  Teuchos::RCP<FlowRelations::EOSVaporInGas> eos_gas_;

  // wrms specified on a region-basis
  typedef std::pair< std::string, Teuchos::RCP<FlowRelations::WRM> > WRMRegionPair;
  std::vector< Teuchos::RCP<WRMRegionPair> > wrm_;

  // ice-liquid capillary pressure model, may be different in ice wedges?
  // for now, just a single one...
  Teuchos::RCP<FlowRelations::PCIceWater> pc_ice_liq_model_;

  // mathematical operators
  Teuchos::RCP<Amanzi::BDFTimeIntegrator> time_stepper_;
  Teuchos::RCP<Operators::MatrixMFD> matrix_;
  Teuchos::RCP<Operators::MatrixMFD> preconditioner_;
  double atol_;
  double rtol_;
  double time_step_reduction_factor_;

  // boundary condition data
  Teuchos::RCP<BoundaryFunction> bc_pressure_;
  Teuchos::RCP<BoundaryFunction> bc_head_;
  Teuchos::RCP<BoundaryFunction> bc_flux_;
  std::vector<Operators::Matrix_bc> bc_markers_;
  std::vector<double> bc_values_;

  // factory registration
  static RegisteredPKFactory<Permafrost> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
