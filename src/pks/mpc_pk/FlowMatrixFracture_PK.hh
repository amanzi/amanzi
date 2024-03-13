/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
           Daniil Svyatskiy
*/

/*!

The process kernel that couples Darcy flow in matrix and fracture network.

Mathematical models
...................
Let subscripts :math:`m` and :math:`f` correspond to matrix and fracture, respectively.
The conceptual PDE model of the stationary coupled matrix-fracture flow is

.. math::
  \begin{array}{l}
  \phi_m \frac{S_{s,m}}{g} \frac{\partial p_m}{\partial t}
  + \boldsymbol{\nabla} \cdot (\rho_l \boldsymbol{q}_m) = Q_m,
  \quad
  \boldsymbol{q}_m = -\frac{\boldsymbol{K}_m}{\mu}
  (\boldsymbol{\nabla} p_m - \rho_l \boldsymbol{g}) \\
  %
  \phi_f \frac{S_{s,f}}{g} \frac{\partial p_f}{\partial t}
  + \boldsymbol{\nabla} \cdot (\rho_l \boldsymbol{q}_f) =
    -\rho_l [[ \tilde{\boldsymbol{q}}_m \cdot \boldsymbol{n} ]],
  \quad
  \boldsymbol{q}_f = -\frac{\boldsymbol{K}_f}{\mu}
  (\boldsymbol{\nabla} p_f - \rho_l \boldsymbol{g}) \\
  %
  \tilde{\boldsymbol{q}}_m \cdot \boldsymbol{n} = \frac{k}{g} (p_f - p_m)
  \end{array}

subject to convential boundary conditions for both matrix and fracture domains expect for
the matrix-fracture boundary where the boundary condition is

.. math::
  \boldsymbol{q}_m \cdot \boldsymbol{n} = \tilde{\boldsymbol{q}}_m \cdot \boldsymbol{n}

Here
:math:`\rho_l` is fluid density [kg/m^3],
:math:`\phi` is porosity [-],
:math:`S_s` is ispecific storage [m],
:math:`p` is aqueous pressure [Pa],
:math:`\boldsymbol{K}` is absolute permeability [m^2] for matrix domain and [m^3] for fracture domain,
:math:`Q_m` is source or sink term,
:math:`\boldsymbol{q}` is the Darcy velocity [m/s] for matrix domain and [m^2/s] for fracture domain,
:math:`k` is effective normal premeability [s^-1],
and
:math:`\boldsymbol{g}` is gravity [:math:`m/s^2`].


Main parameters and sublists
............................
* `"PKs order`" [array(string)] defines user names for two flow PKs. The matrix PK should be
  defined *first*.

* `"time integrator`" [list] defines a generic time integrator used by the cycle driver.

.. code-block:: xml

  <ParameterList name="PKs">  <!-- parent list -->
  <ParameterList name="_COUPLED DARCY FLOW">
    <Parameter name="PKs order" type="Array(string)" value="{_FLOW MATRIX, _FLOW FRACTURE}"/>
    <Parameter name="master PK index" type="int" value="0"/>
    <ParameterList name="time integrator">
      ...
    </ParameterList>
  </ParameterList>
  </ParameterList>

*/

#ifndef AMANZI_DARCY_MATRIX_FRACTURE_PK_HH_
#define AMANZI_DARCY_MATRIX_FRACTURE_PK_HH_

#include "Teuchos_RCP.hpp"

#include "EvaluatorIndependentFunction.hh"
#include "EvaluatorSecondary.hh"
#include "PK_BDF.hh"
#include "PK_MPCStrong.hh"
#include "PK_Factory.hh"
#include "TreeOperator.hh"

namespace Amanzi {

class FlowMatrixFracture_PK : public PK_MPCStrong<PK_BDF> {
 public:
  FlowMatrixFracture_PK(Teuchos::ParameterList& pk_tree,
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

  virtual void FunctionalResidual(double t_old,
                                  double t_new,
                                  Teuchos::RCP<TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> f) override;

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> u, double dt) override;

  // // preconditioner application
  virtual int
  ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) override;

  std::string name() override { return "flow matrix-fracture"; }

  // virtual void CalculateDiagnostics() {};
  Teuchos::RCP<const Teuchos::ParameterList> linear_operator_list_;
  Teuchos::RCP<const Teuchos::ParameterList> preconditioner_list_;
  Teuchos::RCP<Teuchos::ParameterList> ti_list_;

 private:
  const Teuchos::RCP<Teuchos::ParameterList>& glist_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_matrix_, mesh_fracture_;

  Teuchos::RCP<Operators::TreeOperator> op_matrix_, op_preconditioner_;

  // factory registration
  static RegisteredPKFactory<FlowMatrixFracture_PK> reg_;
};

} // namespace Amanzi

#endif
