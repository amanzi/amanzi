/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

The conceptual PDE model for the energy equation is 

.. math::
  \frac{\partial \varepsilon}{\partial t} 
  =
  \boldsymbol{\nabla} \cdot (\kappa \nabla T) -
  \boldsymbol{\nabla} \cdot (\eta_l H_l \boldsymbol{q}_l) + Q

where 
:math:`\varepsilon` is the energy density [:math:`J/m^3`],
:math:`\eta_l` is molar density of liquid [:math:`mol/m^3`],
:math:`Q` is heat source term,
:math:`\boldsymbol{q}_l` is the Darcy velocity [m/s],
:math:`\kappa` is thermal conductivity,
and :math:`H_l` is molar enthalpy of liquid [J/mol].
We define 

.. math::
   \varepsilon = \phi (\eta_l s_l U_l + \eta_g s_g U_g) + 
   (1 - \phi) \rho_r c_r T

where
:math:`s_l` is liquid saturation [-],
:math:`s_g` is gas saturation (water vapor),
:math:`\eta_l` is molar density of liquid [:math:`mol/m^3`],
:math:`\eta_g` is molar density of gas,
:math:`U_l` is molar internal energy of liquid [J/mol],
:math:`U_g` is molar internal energy of gas (water vapor) [J/mol],
:math:`\phi` is porosity [-],
:math:`\rho_r` is rock density [:math:`kg/m^3`],
:math:`c_r` is specific heat of rock [J/kg/K],
and :math:`T` is temperature [K].


Physical models and assumptions
...............................
This list is used to summarize physical models and assumptions, such as
coupling with other PKs.
This list is often generated on a fly by a high-level MPC PK.

* `"vapor diffusion`" [bool] is set up automatically by a high-level PK,
  e.g. by EnergyFlow PK. The default value is `"false`".

* `"eos lookup table`" [string] provides the name for optional EOS lookup table.

.. code-block:: xml

  <ParameterList>  <!-- parent list -->
  <ParameterList name="_ENERGY"> 
    <ParameterList name="physical models and assumptions">
      <Parameter name="vapor diffusion" type="bool" value="false"/>
      <Parameter name="eos lookup table" type="string" value="h2o.eos"/>
    </ParameterList>
  </ParameterList>
  </ParameterList>

*/

#ifndef AMANZI_ENERGY_PK_HH_
#define AMANZI_ENERGY_PK_HH_

// TPLs
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "EvaluatorPrimary.hh"
#include "Key.hh"
#include "Operator.hh"
#include "PDE_Accumulation.hh"
#include "PDE_AdvectionUpwind.hh"
#include "PDE_Diffusion.hh"
#include "PK.hh"
#include "PK_BDF.hh"
#include "PK_DomainFunction.hh"
#include "PK_PhysicalBDF.hh"
#include "Tensor.hh"
#include "TreeVector.hh"
#include "Upwind.hh"
#include "VerboseObject.hh"

namespace Amanzi {
namespace Energy {

class Energy_PK : public PK_PhysicalBDF {
 public:
  Energy_PK(Teuchos::ParameterList& pk_tree,
            const Teuchos::RCP<Teuchos::ParameterList>& glist,
            const Teuchos::RCP<State>& S,
            const Teuchos::RCP<TreeVector>& soln);
  virtual ~Energy_PK(){};

  // methods required by PK interface
  virtual void Setup() override;
  virtual void Initialize() override;
  virtual std::string name() override { return passwd_; }

  // methods required for time integration
  // -- management of the preconditioner
  virtual int
  ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> hu) override
  {
    return op_preconditioner_->ApplyInverse(*u->Data(), *hu->Data());
  }

  // -- check the admissibility of a solution
  //    override with the actual admissibility check
  bool IsAdmissible(Teuchos::RCP<const TreeVector> up) override { return true; }

  // -- possibly modifies the predictor that is going to be used as a
  //    starting value for the nonlinear solve in the time integrator,
  //    the time integrator will pass the predictor that is computed
  //    using extrapolation and the time step that is used to compute
  //    this predictor this function returns true if the predictor was
  //    modified, false if not
  bool
  ModifyPredictor(double dt, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u) override
  {
    return false;
  }

  // -- possibly modifies the correction, after the nonlinear solver (NKA)
  //    has computed it, will return true if it did change the correction,
  //    so that the nonlinear iteration can store the modified correction
  //    and pass it to NKA so that the NKA space can be updated
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double dt,
                   Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<TreeVector> du) override;

  // -- calling this indicates that the time integration
  //    scheme is changing the value of the solution in state.
  void ChangedSolution() override { temperature_eval_->SetChanged(); }

  // other methods
  bool UpdateConductivityData(const Teuchos::Ptr<State>& S);
  void UpdateSourceBoundaryData(double T0, double T1, const CompositeVector& u);
  void ComputeBCs(const CompositeVector& u);
  void AddSourceTerms(CompositeVector& rhs);

  // access
  virtual Teuchos::RCP<Operators::Operator>
  my_operator(const Operators::OperatorType& type) override;

  virtual Teuchos::RCP<Operators::PDE_HelperDiscretization>
  my_pde(const Operators::PDEType& type) override
  {
    return op_matrix_diff_;
  }

  // -- for unit tests
  std::vector<WhetStone::Tensor>& get_K() { return K; }

 protected:
  void InitializeFieldFromField_(const std::string& field0,
                                 const std::string& field1,
                                 bool call_evaluator);

 private:
  void InitializeFields_();

 public:
  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;

 protected:
  int dim;

  const Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::RCP<Teuchos::ParameterList> ep_list_;
  Teuchos::RCP<const Teuchos::ParameterList> preconditioner_list_;
  Teuchos::RCP<Teuchos::ParameterList> ti_list_;

  // primary field
  std::string passwd_;
  Teuchos::RCP<EvaluatorPrimary<CompositeVector, CompositeVectorSpace>> temperature_eval_;

  // names of state fields
  Key temperature_key_;
  Key energy_key_, prev_energy_key_;
  Key enthalpy_key_, aperture_key_, prev_aperture_key_;
  Key ie_liquid_key_, ie_gas_key_, ie_rock_key_;
  Key vol_flowrate_key_, particle_density_key_, sat_liquid_key_;
  Key mol_density_liquid_key_, mass_density_liquid_key_;
  Key mol_density_gas_key_, x_gas_key_;
  Key conductivity_gen_key_, conductivity_key_, conductivity_eff_key_;

  // conductivity tensor
  std::vector<WhetStone::Tensor> K;

  // boundary conditons
  std::vector<Teuchos::RCP<PK_DomainFunction>> bc_temperature_;
  std::vector<Teuchos::RCP<PK_DomainFunction>> bc_flux_;
  int dirichlet_bc_faces_;

  // source terms
  std::vector<Teuchos::RCP<PK_DomainFunction>> srcs_;

  // operators and solvers
  Teuchos::RCP<Operators::PDE_Diffusion> op_matrix_diff_, op_preconditioner_diff_;
  Teuchos::RCP<Operators::PDE_Accumulation> op_acc_;
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> op_matrix_advection_, op_preconditioner_advection_;
  Teuchos::RCP<Operators::Operator> op_matrix_, op_preconditioner_, op_advection_;
  Teuchos::RCP<Operators::BCs> op_bc_, op_bc_enth_;

  bool prec_include_enthalpy_;

  // upwinding
  Teuchos::RCP<CompositeVector> upw_conductivity_;
  Teuchos::RCP<Operators::Upwind> upwind_; // int implies fake model

  // fracture network
  bool flow_on_manifold_;
};

} // namespace Energy
} // namespace Amanzi

#endif
