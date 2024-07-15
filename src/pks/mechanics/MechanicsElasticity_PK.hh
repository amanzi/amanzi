/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

The conceptual PDE model for the quasi-static elastic deformation is

.. math::
  \boldsymbol{\nabla} \cdot (C \varepsilon(\boldsymbol{d}))
  =
  \rho \boldsymbol{g}

where 
:math:`\boldsymbol{d}` is the displacement [m],
:math:`\rho` is the rock density [kg/m^3],
:math:`\boldsymbol{\varepsilon}` is the strain tensor,
and
:math:`\boldsymbol{g}` is the gravity vector [:math:`m/s^2`].

For a linear elasticity problem, the stress tensor :math:`C` is a linear operator 
acting on the strain tensor.
In general, the stress tensor is a nonlinear operator. 


Global parameters
.................
Global parameters are placed in the sublist `"mechanics`".
The list of global parameters include:

.. admonition:: mechanics_pk-spec

  * `"domain name`" ``[string]`` specifies mesh name that defined domain of this PK.
    Default is `"domain`".


Physical models and assumptions
...............................
This list is used to summarize physical models and assumptions, such as
coupling with other PKs.
This list is often generated or extended by a high-level MPC PK.

.. admonition:: mechanics_assumptions-spec

  * `"use gravity`" ``[bool]`` defines non-zero source term. Default is *false*.

  * `"biot scheme: undrained split`" ``[bool]`` defines iterative coupling with 
    a flow PK where mechanics is solved first.

  * `"biot scheme: fixed-stress split`" ``[bool]`` defines iterative coupling with 
    a flow PK where flow is solved first.

.. code-block:: xml

  <ParameterList name="mechanics">  <!-- parent list -->
  <ParameterList name="physical models and assumptions">
    <Parameter name="use gravity" type="bool" value="false"/>
    <Parameter name="biot scheme: undrained split" type="bool" value="false"/>
    <Parameter name="biot scheme: fixed stress split" type="bool" value="false"/>
  </ParameterList>


Main sublists
.............
The following sublists are needed to create and control this PK:

.. code-block:: xml

  <ParameterList name="mechanics">  <!-- parent list -->
    <ParameterList name="time integrator">
    <ParameterList name="operators">
    <ParameterList name="physical models and assumptions">
    <ParameterList name="boundary conditions">

*/

#ifndef AMANZI_MECHANICS_ELASTICITY_PK_HH_
#define AMANZI_MECHANICS_ELASTICITY_PK_HH_

#include <string>
#include <vector>

// TPLs
#include "Epetra_Vector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "BDF1_TI.hh"
#include "EvaluatorPrimary.hh"
#include "Key.hh"
#include "PDE_Accumulation.hh"
#include "PDE_Abstract.hh"
#include "PDE_Elasticity.hh"
#include "PK.hh"
#include "PK_Factory.hh"
#include "PK_PhysicalBDF.hh"
#include "State.hh"
#include "TreeOperator.hh"
#include "TreeVector.hh"
#include "Units.hh"
#include "VerboseObject.hh"

#include "HydrostaticStressEvaluator.hh"
#include "Mechanics_PK.hh"
#include "MechanicsBoundaryFunction.hh"
#include "VolumetricStrainEvaluator.hh"

namespace Amanzi {
namespace Mechanics {

class MechanicsElasticity_PK : public Mechanics_PK {
 public:
  MechanicsElasticity_PK(Teuchos::ParameterList& pk_tree,
                         const Teuchos::RCP<Teuchos::ParameterList>& glist,
                         const Teuchos::RCP<State>& S,
                         const Teuchos::RCP<TreeVector>& soln);

  ~MechanicsElasticity_PK(){};

  // methods required for PK interface
  virtual void Setup() final;
  virtual void Initialize() final;

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) final;
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) final;

  virtual std::string name() override { return "mechanics"; }

  // methods required for time integration interface
  // -- computes the non-linear functional f = f(t,u,udot) and related norm.
  void FunctionalResidual(const double t_old,
                          double t_new,
                          Teuchos::RCP<TreeVector> u_old,
                          Teuchos::RCP<TreeVector> u_new,
                          Teuchos::RCP<TreeVector> f) override;
  double ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) override;

  // -- management of the preconditioner
  int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> pu) override;
  void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> u, double dt) override;

  // -- access
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace>> bdf1_dae() { return bdf1_dae_; }
  virtual Teuchos::RCP<Operators::Operator>
  my_operator(const Operators::OperatorType& type) override;
  virtual Teuchos::RCP<Operators::PDE_HelperDiscretization>
  my_pde(const Operators::PDEType& type) override
  {
    return op_matrix_elas_;
  }

 public:
  Teuchos::RCP<Operators::Operator> op_matrix_;
  Teuchos::RCP<Operators::PDE_Elasticity> op_matrix_elas_, op_matrix_graddiv_;

 private:
  static RegisteredPKFactory<MechanicsElasticity_PK> reg_;
};

} // namespace Mechanics
} // namespace Amanzi

#endif
