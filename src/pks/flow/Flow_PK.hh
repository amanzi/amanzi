/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (version 1)
           Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

/*!

Mathematical models
...................
A few PDE models can be instantiated using the parameters described below.


Fully saturated flow
````````````````````
The conceptual PDE model for the fully saturated flow is

.. math::
  \left(\frac{S_s}{g} + \frac{S_y}{Lg}\right)\frac{\partial p_l}{\partial t}
  =
  -\boldsymbol{\nabla} \cdot (\rho_l \boldsymbol{q}_l) + Q,
  \quad
  \boldsymbol{q}_l
  = -\frac{\boldsymbol{K}}{\mu}
  (\boldsymbol{\nabla} p - \rho_l \boldsymbol{g}),

where
:math:`S_s` and :math:`S_y` are specific storage [1/m] and specific yield [-], respectively,
:math:`L` is characteristic length [m],
:math:`\rho_l` is fluid density [:math:`kg / m^3`],
:math:`Q` is source or sink term [:math:`kg / m^3 / s`],
:math:`\boldsymbol{q}_l` is the Darcy velocity [:math:`m/s`],
and :math:`\boldsymbol{g}` is gravity [:math:`m/s^2`].
The specific storgae can be defined using

.. math::
  S_s = \left(\phi\, \beta_f + \beta_m\right)\rho_l\,g

where :math:`\beta_f` [1/Pa] and :math:`\beta_m` [1/Pa] are fluid and matrix compressibilities, respectively.


Partially saturated flow with water vapor
`````````````````````````````````````````
The conceptual PDE model for the partially saturated flow with water vapor
includes liquid phase (liquid water) and gas phase (water vapor):

.. math::
  \frac{\partial \theta}{\partial t}
  =
  - \boldsymbol{\nabla} \cdot (\eta_l \boldsymbol{q}_l)
  - \boldsymbol{\nabla} \cdot (\boldsymbol{K}_g \boldsymbol{\nabla} \big(\frac{p_v}{p_g}\big)) + Q,
  \quad
  \boldsymbol{q}_l
  = -\frac{\boldsymbol{K} k_r}{\mu}
  (\boldsymbol{\nabla} p - \rho_l \boldsymbol{g})

where
:math:`\theta` is total water storage [:math:`mol/m^3`],
:math:`\eta_l` is molar density of liquid (water) [:math:`mol/m^3`],
:math:`\rho_l` is fluid density [:math:`kg/m^3`],
:math:`Q` is source or sink term [:math:`mol/m^3/s`],
:math:`\boldsymbol{q}_l` is the Darcy velocity [:math:`m/s`],
:math:`k_r` is relative permeability [-],
:math:`\boldsymbol{g}` is gravity [:math:`m/s^2`],
:math:`p_v` is the vapor pressure [Pa],
:math:`p_g` is the gas pressure [Pa],
and :math:`\boldsymbol{K}_g` is the effective diffusion coefficient of the water vapor.
We define

.. math::
  \theta = \phi \eta_l s_l + \phi \eta_g (1 - s_l) X_g

where :math:`s_l` is liquid saturation [-],
:math:`\phi` is porosity [-],
:math:`\eta_g` is molar density of water vapor [:math:`mol/m^3`],
and :math:`X_g` is molar fraction of water vapor.
The effective diffusion coefficient of the water vapor is given by

.. math::
  \boldsymbol{K}_g = \phi s_g \tau_g \eta_g \boldsymbol{D}_g

where :math:`s_g` is gas saturation [-],
:math:`\tau_g` is the tortuosity of the gas phase [-],
:math:`\eta_g` is the molar density of gas [:math:`kg/m^3`],
and :math:`\boldsymbol{D}_g` is the diffusion coefficient of the gas phase [:math:`m^2/s`],
The gas pressure :math:`p_g` is set to the atmosperic pressure and the vapor pressure
model assumes thermal equlibrium of liquid and gas phases:

.. math::
  p_v = P_{sat}(T) \exp\left(\frac{P_{cgl}}{\eta_l R T}\right)

where
:math:`R` is the ideal gas constant [:math:`kg m^2/K/mol/s^2`],
:math:`P_{cgl}` is the liquid-gas capillary pressure [Pa],
:math:`P_{sat}` is the saturated vapor pressure [Pa],
and :math:`T` is the temperature [K].
The diffusion coefficient is based of TOUGHT2 model

.. math::
   D_g = D_0 \frac{P_{ref}}{p} \left(\frac{T}{273.15}\right)^a

where
:math:`D_0 = 2.14 \cdot 10^{-5}`,
:math:`P_{ref}` is atmospheric pressure,
and :math:`a = 1.8`.
finally we need a model for the gas tortuosity. We use the Millington and Quirk model:

.. math::
   \tau_g = \phi^\beta s_g^\gamma

where
:math:`\beta = 1/3` and
:math:`\gamma = 7/3`.


Isothermal partially saturated flow with dual porosity model
````````````````````````````````````````````````````````````
The conceptual model for the partially saturated flow with dual porosity model
assumes that water flow is restricted to the fractures and the water in the matrix does not move.
The rock matrix represents immobile pockets that can exchange, retain and store water
but do not permit convective flow.
This leads to dual-porosity type flow and transport models that partition the liquid
phase into mobile and immobile regions.
The Richards equation in the mobile region is augmented by the water exchange
term :math:`\Sigma_w`:

.. math::
  \frac{\partial \theta_{lf}}{\partial t}
  = -\boldsymbol{\nabla} \cdot (\eta_l \boldsymbol{q}_l)
    -\frac{K_m\,k_{rm}\,\eta_l}{\mu\,L_m}\, \nabla p_m + Q_f,
  \qquad
  \boldsymbol{q}_l
  = -\frac{\boldsymbol{K}_f\, k_{rf}}{\mu}
  (\boldsymbol{\nabla} p_f - \rho_l \boldsymbol{g})

where
:math:`p_f` is fracture pressure [Pa],
:math:`p_m` is matrix pressure [Pa],
:math:`L_m` is the characteristic matrix depth defined typically as the ratio of a matrix block [m],
and :math:`Q_f` is source or sink term [:math:`kg \cdot m^{-3} \cdot s^{-1}`].
The equation for water balance in the matrix is

.. math::
  \frac{\partial \theta_{lm}}{\partial t}
  = Q_m
    +\nabla\cdot \left(\frac{K_m\, k_{rm}\,\eta_l}{\mu}\, \nabla p_{m}\right),

where
:math:`Q_m` is source or sink term [:math:`kg / m^3 / s`].
The water storages are defined as

.. math::
  \theta_f = \phi_f\, \eta_l\, s_{lf},\quad
  \theta_m = \phi_m\, \eta_l\, s_{lm},

where saturations :math:`s_{lf}` and :math:`s_{lm}` may use different capillary
pressure - saturation models.
In the simplified model, the rate of water exchange between the fracture and matrix regions
is proportional to the difference in hydraulic heads:

.. math::
  \frac{K_m\,k_{rm}\,\eta_l}{\mu\,L_m}\, \nabla p_m
  \approx
  \alpha_w (h_f - h_m),

where :math:`\alpha_w` is the mass transfer coefficient.
Since hydraulic heads are needed for both regions, this equation requires
retention curves for both regions and therefore is nonlinear.


Physical models and assumptions
...............................
This list is used to summarize physical models and assumptions, such as
coupling with other PKs.
This list is often generated or extended by a high-level MPC PK.
In the code development, this list plays a two-fold role.
First, it provides necessary information for coupling different PKs such
as flags for adding a vapor diffusion to Richards' equation.
Second, the developers may use it instead of a factory of evaluators such as
creation of primary and secondary evaluators for rock porosity models.
Combination of both approaches may lead to a more efficient code.

* `"vapor diffusion`" [bool] is set up automatically by a high-level PK,
  e.g. by EnergyFlow PK. The default value is `"false`".

* `"flow and transport in fractures`" [bool] indicates that Darcy flow is calculated in fractures.
  This option is ignored is mesh dimentionaly equals to manifold dimensionality.

* `"multiscale model`" [string] specifies a multiscale model.
  Available options are `"single porosity`" (default) and `"dual continuum discontinum matrix`".

* `"viscosity model`" [string] changes the evaluator for liquid viscosity.
  Available options are `"generic`" and `"constant viscosity`" (default).

* `"porosity model`" [string] specifies an isothermal porosity model.
  Available options are `"compressible: storativity coefficient`",
  `"compressible: pressure function`", and `"constant porosity`" (default).

* `"coupled matrix fracture flow`" [string] specifies PK's role in the strong
  coupling of two flow PKs. The value is either `"matrix`" or `"fracture`".

* `"eos lookup table`" [string] provides the name for optional EOS lookup table.

.. code-block:: xml

  <ParameterList name="flow">  <!-- parent list -->
  <ParameterList name="physical models and assumptions">
    <Parameter name="vapor diffusion" type="bool" value="false"/>
    <Parameter name="viscosity model" type="string" value="constant viscosity"/>
    <Parameter name="porosity model" type="string" value="compressible: pressure function"/>
    <Parameter name="multiscale model" type="string" value="single porosity"/>
    <Parameter name="coupled matrix fracture flow" type="string" value="matrix"/>
    <Parameter name="eos lookup table" type="string" value="h2o.eos"/>
  </ParameterList>
  </ParameterList>


Global parameters
.................

* `"domain name`" [string] specifies mesh name that defined domain of this PK.
  Default is `"domain`".

*/

#ifndef AMANZI_FLOW_PK_HH_
#define AMANZI_FLOW_PK_HH_

#include <vector>

// TPLs
#include "Epetra_IntVector.h"
#include "Epetra_FECrsMatrix.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "BCs.hh"
#include "BDFFnBase.hh"
#include "Checkpoint.hh"
#include "CompositeVectorSpace.hh"
#include "EvaluatorIndependentFunction.hh"
#include "Key.hh"
#include "Operator.hh"
#include "PK_DomainFunction.hh"
#include "PK_PhysicalBDF.hh"
#include "EvaluatorPrimary.hh"
#include "Tensor.hh"
#include "Units.hh"
#include "VerboseObject.hh"

// Flow
#include "FlowBoundaryFunction.hh"
#include "FlowSourceFunction.hh"
#include "FlowDefs.hh"
#include "FlowTypeDefs.hh"

namespace Amanzi {
namespace Flow {

class Flow_PK : public PK_PhysicalBDF {
 public:
  Flow_PK();
  Flow_PK(Teuchos::ParameterList& pk_tree,
          const Teuchos::RCP<Teuchos::ParameterList>& glist,
          const Teuchos::RCP<State>& S,
          const Teuchos::RCP<TreeVector>& soln);
  virtual ~Flow_PK(){};

  // members required by PK interface
  virtual void Setup() override;
  virtual void Initialize() override;

  // other members of this PK.
  // -- initialize simple fields common for both flow models.
  void UpdateLocalFields_(const Teuchos::Ptr<State>& S);

  // --- management of boundary and source terms
  void UpdateSourceBoundaryData(double t_old, double t_new, const CompositeVector& u);
  void ComputeOperatorBCs(const CompositeVector& u);
  void SeepageFacePFloTran(const CompositeVector& u, int* nseepage, double* area_seepage);
  void SeepageFaceFACT(const CompositeVector& u, int* nseepage, double* area_seepage);

  void AddSourceTerms(CompositeVector& rhs, double dt);
  void ComputeWellIndex(Teuchos::ParameterList& spec);
  bool IsWellIndexRequire(Teuchos::ParameterList& spec);

  // -- absolute permeability and derived quantities.
  void SetAbsolutePermeabilityTensor();

  // -- miscallenous members
  void DeriveFaceValuesFromCellValues(const Epetra_MultiVector& ucells, Epetra_MultiVector& ufaces);

  // -- io members
  void OutputTimeHistory(const Teuchos::ParameterList& plist, std::vector<dt_tuple>& dt_history);

  // -- utilities
  double WaterVolumeChangePerSecond(const std::vector<int>& bc_model,
                                    const Epetra_MultiVector& vol_flowrate) const;

  // -- V&V
  void VV_ValidateBCs() const;
  void VV_ReportWaterBalance(const Teuchos::Ptr<State>& S) const;
  void VV_ReportSeepageOutflow(const Teuchos::Ptr<State>& S, double dT) const;
  void VV_PrintHeadExtrema(const CompositeVector& pressure) const;
  void VV_PrintSourceExtrema() const;
  void VV_FractureConservationLaw() const;

  // -- extensions
  void VerticalNormals(int c, AmanziGeometry::Point& n1, AmanziGeometry::Point& n2);
  virtual double BoundaryFaceValue(int f, const CompositeVector& u);

  // access
  Teuchos::RCP<Operators::BCs> op_bc() { return op_bc_; }
  double seepage_mass() { return seepage_mass_; } // support of unit tests

 private:
  void InitializeFields_();

 protected:
  void InitializeBCsSources_(Teuchos::ParameterList& list);

 public:
  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;

  double dt_, dt_next_;

  int MyPID; // parallel information: will be moved to private
  int missed_bc_faces_, dirichlet_bc_faces_;
  int ti_phase_counter;

 public:
  Teuchos::RCP<const Teuchos::ParameterList> linear_operator_list_;
  Teuchos::RCP<const Teuchos::ParameterList> preconditioner_list_;
  Teuchos::RCP<Teuchos::ParameterList> ti_list_;

 protected:
  int dim;

  std::string passwd_;
  bool peaceman_model_, use_bulk_modulus_;

  // Stationary physical quantatities
  std::vector<WhetStone::Tensor> K;
  AmanziGeometry::Point gravity_;
  double g_, rho_, molar_rho_, atm_pressure_;
  double flux_units_; // scaling for flux units from kg to moles.

  Teuchos::RCP<Epetra_MultiVector> Kxy;
  std::string coordinate_system_;

  // boundary conditions
  std::vector<Teuchos::RCP<FlowBoundaryFunction>> bcs_;
  int nseepage_prev;

  Teuchos::RCP<Operators::BCs> op_bc_;

  // source terms and liquid balance
  std::vector<Teuchos::RCP<FlowSourceFunction>> srcs;
  mutable double mass_bc, seepage_mass_, mass_initial;

  // field evaluators (MUST GO AWAY lipnikov@lanl.gov)
  Teuchos::RCP<EvaluatorPrimary<CompositeVector, CompositeVectorSpace>> vol_flowrate_eval_;
  Teuchos::RCP<EvaluatorPrimary<CompositeVector, CompositeVectorSpace>> pressure_eval_,
    pressure_msp_eval_;

  // DFN model
  bool flow_on_manifold_; // true for the DFN model
  bool coupled_to_matrix_, coupled_to_fracture_;

  // names of state fields
  Key pressure_key_;
  Key vol_flowrate_key_, mol_flowrate_key_, darcy_velocity_key_;
  Key wc_key_, specific_storage_key_, specific_yield_key_;
  Key saturation_liquid_key_, prev_saturation_liquid_key_;
  Key porosity_key_, hydraulic_head_key_, pressure_head_key_;
  Key permeability_key_, permeability_eff_key_;
  Key water_storage_key_, prev_water_storage_key_;
  Key viscosity_liquid_key_, mol_density_liquid_key_;
  Key prev_aperture_key_, aperture_key_, bulk_modulus_key_;

  // io
  Utils::Units units_;
  Teuchos::RCP<Teuchos::ParameterList> fp_list_;
};

} // namespace Flow
} // namespace Amanzi

#endif
