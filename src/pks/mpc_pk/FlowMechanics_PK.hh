/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*!

Sequential coupling of fully flow and mechanics PKs via the fixed stress
split algorithm.
The conceptual PDE model of the coupled flow and mechanics equations is

.. math::
  \begin{array}{l}
  -\boldsymbol{\nabla} \cdot (C \varepsilon(\boldsymbol{d})) + \alpha \nabla p + \beta_r \nabla T
  =
  \rho \boldsymbol{g} \\[1ex]
  %
  \displaystyle\frac{\partial (\eta \phi)}{\partial t}
  =
  - \boldsymbol{\nabla} \cdot (\eta \boldsymbol{q}) + Q_1\\[1ex]
  %
  \displaystyle\frac{\partial (\phi\, \eta\, U_l + (1 - \phi) \rho_r U_r)}{\partial t}
  =
  \boldsymbol{\nabla} \cdot (\kappa \nabla T) -
  \boldsymbol{\nabla} \cdot (\eta\, H\, \boldsymbol{q}) + Q_2.
  \end{array}

In the first equation,
:math:`C` is elasticity tensor [Pa],
:math:`\varepsilon` is stain tensor [-],
:math:`\boldsymbol{d}` is displacement [m],
:math:`\rho` is density of liquid [:math:`kg/m^3`],
:math:`\alpha` is Biot coefficient [-],
:math:`\beta_r` is thermal stress coefficient [:math:`Pa/K`],
:math:`\boldsymbol{g}` is gravity [:math:`m/s^2`],
:math:`p` is pressure [Pa],
and
:math:`T` is temperature [K].

In the second equation,
:math:`\eta` is molar density of liquid [:math:`mol/m^3`],
:math:`\phi` is porosity [-],
:math:`Q_1` is source or sink term,
and
:math:`\boldsymbol{q}` is Darcy velocity [m/s].
The porosity model is

.. math::
  \phi = \phi_0 + \alpha\, {\rm div}\, \varepsilon + c_0\,(p - p_0) - \beta_l (T - T_0).

where
:math:`\beta_l` is thermal dilation coefficient [:math:`Pa/K`],
:math:`c_0` is pore compressibility [:math:`Pa^{-1}`],
:math:`p_0` is reference pressure [:math:`Pa`],
and
:math:`T_0` is reference temperature [:math:`K`].
Currently :math:`T_0 = 0^\circ C`.

In the third equation,
:math:`\kappa` is thermal conductivity,
:math:`U_l` is internal energy of liquid [J/mol],
:math:`U_r` is internal energy of rock [J/kg],
:math:`H` is molar enthalpy of liquid [J/mol],
:math:`\rho_r` is rock density [:math:`kg/m^3`],
and
:math:`Q_2` is heat source term.
The internal energies are based on linear models:

.. math::
  U_i = c_i\,(T - T_0),\quad i=l,r,

where :math:`c_i` are specific heat capacities.

*/

#ifndef AMANZI_FLOW_MECHANICS_PK_HH_
#define AMANZI_FLOW_MECHANICS_PK_HH_

#include "Key.hh"

#include "PK_MPCSequential.hh"

namespace Amanzi {

class FlowMechanics_PK : public PK_MPCSequential {
 public:
  FlowMechanics_PK(Teuchos::ParameterList& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  virtual void Setup() override;
  virtual void Initialize() override;
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) override;
  virtual double
  ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) override;
  virtual void CommitSequentialStep(Teuchos::RCP<const TreeVector> u_old,
                                    Teuchos::RCP<const TreeVector> u_new) override;

 private:
  Key domain_;
  const Teuchos::RCP<Teuchos::ParameterList> glist_;

  Key displacement_key_, hydrostatic_stress_key_, vol_strain_key_, biot_key_;
  Key pressure_key_, porosity_key_, saturation_liquid_key_, water_storage_key_;

  bool thermal_flow_;

 private:
  static RegisteredPKFactory<FlowMechanics_PK> reg_;
};

} // namespace Amanzi

#endif
