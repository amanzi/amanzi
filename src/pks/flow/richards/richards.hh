/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
A base two-phase, thermal Richard's equation with water vapor.

License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#ifndef PK_FLOW_RICHARDS_HH_
#define PK_FLOW_RICHARDS_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "composite_vector.hh"
#include "tree_vector.hh"
#include "state.hh"
#include "matrix_mfd.hh"

#include "PK.hh"
#include "bdf_fn_base.hh"
#include "bdf_time_integrator.hh"
#include "flow-bc-factory.hh"

#include "water_retention_model.hh"
#include "eos.hh"
#include "flow.hh"

namespace Amanzi {
namespace Flow {

class Richards : public Flow {

public:
  Richards(Teuchos::ParameterList& flow_plist, const Teuchos::RCP<State>& S,
           const Teuchos::RCP<TreeVector>& solution);

  // main methods
  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::RCP<State>& S);

  // -- Advance from state S to state S_next at time S0.time + dt.
  virtual bool advance(double dt);

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) = 0;

  // applies preconditioner to u and returns the result in Pu
  virtual void precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);
  
private:

  // computational/convenience methods
  void DeriveDarcyFlux_(const Epetra_Vector& flux, Epetra_MultiVector& velocity);
  void DeriveDarcyVelocity_(const Epetra_Vector& flux, Epetra_MultiVector& velocity);

  void SetAbsolutePermeabilityTensor_(const Teuchos::RCP<State>& S,
          std::vector<WhetStone::Tensor>& K);

  void CalculateRelativePermeabilityUpwindGravity_(const Teuchos::RCP<State>& S,
          const Epetra_Vector& pres, Teuchos::RCP<CompositeVector>& rel_perm);
  void CalculateRelativePermeabilityUpwindFlux_(const Epetra_Vector& pres,
          const Epetra_Vector& flux, Teuchos::RCP<CompositeVector>& rel_perm);
  void CalculateRelativePermeabilityArithmeticMean_(const Teuchos::RCP<State>& S,
          const Epetra_Vector& pres, Teuchos::RCP<CompositeVector>& rel_perm);

  // physical methods
  void ApplyDiffusion_(const Teuchos::RCP<State>& S,const Teuchos::RCP<CompositeVector>& g);
  void AddAccumulation_(const Teuchos::RCP<CompositeVector>& g);
  void UpdateSecondaryVariables_(const Teuchos::RCP<State>& S);

  void DensityLiquid_(const Teuchos::RCP<State>& S, const CompositeVector& temp,
                      const CompositeVector& pres,
                      const Teuchos::RCP<CompositeVector>& dens_liq,
                      const Teuchos::RCP<CompositeVector>& mol_dens_liq);

  void DensityGas_(const Teuchos::RCP<State>& S, const CompositeVector& temp,
                   const CompositeVector& pres, const double& p_atm,
                   const Teuchos::RCP<CompositeVector>& mol_frac_gas,
                   const Teuchos::RCP<CompositeVector>& dens_gas,
                   const Teuchos::RCP<CompositeVector>& mol_dens_gas);

  void ViscosityLiquid_(const Teuchos::RCP<State>& S, const CompositeVector& temp,
                        const Teuchos::RCP<CompositeVector>& visc_liq);
  
  void Saturation_(const Teuchos::RCP<State>& S, const CompositeVector& pres,
                   const double& p_atm, const Teuchos::RCP<CompositeVector>& sat_liq);

  void RelativePermeability_(const Teuchos::RCP<State>& S,
                             const CompositeVector& pres, const double& p_atm,
                             const Teuchos::RCP<CompositeVector>& rel_perm);
  
 private:
  Teuchos::RCP<Amanzi::BDFTimeIntegrator> time_stepper_;
  int Krel_method_;

  std::vector<Teuchos::RCP<AmanziFlow::WaterRetentionModel> > wrms_;
  Teuchos::RCP<FlowRelations::EOS> eos_liquid_;
  Teuchos::RCP<FlowRelations::EOS> eos_gas_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

