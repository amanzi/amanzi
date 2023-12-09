/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
           Konstantin Lipnikov
*/

/*!

Process kernel that strongly couples Flow PK with Energy PK.
The conceptual PDE model of the coupled flow and energy equations is

.. math::
  \begin{array}{l}
  \frac{\partial \theta}{\partial t}
  =
  - \boldsymbol{\nabla} \cdot (\eta_l \boldsymbol{q}_l)
  - \boldsymbol{\nabla} \cdot (\phi s_g \tau_g D_g \boldsymbol{\nabla} X_g) + Q_1,
  \quad
  \boldsymbol{q}_l
  = -\frac{\boldsymbol{K} k_r}{\mu}
  (\boldsymbol{\nabla} p - \rho_l \boldsymbol{g}) \\
  %
  \frac{\partial \varepsilon}{\partial t}
  =
  \boldsymbol{\nabla} \cdot (\kappa \nabla T) -
  \boldsymbol{\nabla} \cdot (\eta_l H_l \boldsymbol{q}_l) + Q_2
  \end{array}

In the first equation,
:math:`\theta` is total water storage (we use non-conventional definition) [:math:`mol/m^3`],
:math:`\eta_l` is molar density of liquid [:math:`mol/m^3`],
:math:`\rho_l` is fluid density [:math:`kg/m^3`],
:math:`Q_1` is source or sink term,
:math:`\boldsymbol{q}_l` is the Darcy velocity [m/s],
:math:`k_r` is relative permeability [-],
:math:`\boldsymbol{g}` is gravity [:math:`m/s^2`],
:math:`\phi` is porosity [-],
:math:`s_g` is gas saturation (water vapor) [-],
:math:`\tau_g` is tortuosity of gas [-],
:math:`D_g` is diffusion coefficient,
and :math:`X_g` is molar fraction of water in the gas phase [-].
We define

.. math::
   \theta = \phi (s_g \eta_g X_g + s_l \eta_l)

where
:math:`s_l` is liquid saturation [-],
and :math:`\eta_g` is molar density of gas.

In the second equation,
:math:`\varepsilon` is the energy density [:math:`J/mol^3`],
:math:`Q_2` is source or sink term,
:math:`\kappa` is thermal conductivity [W/m/K],
:math:`H_l` is molar enthalphy of liquid [J/mol],
and :math:`T` is temperature [K].
We define

.. math::
   \varepsilon = \phi (\eta_l s_l U_l + \eta_g s_g U_g) +
   (1 - \phi) \rho_r c_r T

where
:math:`U_l` is molar internal energy of liquid [J/mol],
:math:`U_g` is molar internal energy of gas (water vapor) [J/mol],
:math:`\rho_r` is rock density [kg/m^3],
and :math:`c_r` is specific heat of rock [J/kg/K].

*/

#ifndef AMANZI_FLOW_ENERGY_PK_HH_
#define AMANZI_FLOW_ENERGY_PK_HH_

#include "Teuchos_RCP.hpp"

#include "EvaluatorIndependentFunction.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "PK_BDF.hh"
#include "PK_MPCStrong.hh"
#include "PK_Factory.hh"

namespace Amanzi {

class FlowEnergy_PK : public PK_MPCStrong<PK_BDF> {
 public:
  FlowEnergy_PK(Teuchos::ParameterList& pk_tree,
                const Teuchos::RCP<Teuchos::ParameterList>& glist,
                const Teuchos::RCP<State>& S,
                const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  virtual void Setup() override;
  virtual void Initialize() override;

  // -- dt is the minimum of the sub pks
  // virtual double get_dt();
  // virtual void set_dt(double dt);

  // -- advance each sub pk from t_old to t_new.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;

  // -- computes non-linear functional f = f(t,u)
  virtual void FunctionalResidual(double t_old,
                                  double t_new,
                                  Teuchos::RCP<TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> f) override;

  // -- enorm for the coupled system
  std::string name() override { return "thermal flow"; }

 private:
  const Teuchos::RCP<Teuchos::ParameterList>& glist_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Key domain_; // computational domain

  Teuchos::RCP<EvaluatorIndependentFunction> particle_density_eval;
  Teuchos::RCP<EvaluatorIndependentFunction> porosity_eval;
  Teuchos::RCP<EvaluatorIndependentFunction> saturation_liquid_eval;

  // keys
  Key pressure_key_, temperature_key_;
  Key ie_rock_key_, ie_gas_key_, ie_liquid_key_, energy_key_, prev_energy_key_;
  Key particle_density_key_;
  Key mol_density_liquid_key_, mol_density_gas_key_, mass_density_liquid_key_;
  Key sat_liquid_key_, prev_sat_liquid_key_, viscosity_liquid_key_;
  Key wc_key_, prev_wc_key_;

  // eos
  std::string eos_table_;

  // factory registration
  static RegisteredPKFactory<FlowEnergy_PK> reg_;
};

} // namespace Amanzi
#endif
