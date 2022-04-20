/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Evaluator_OperatorApply is a secondary evaluator that calculates r = b - Ax

/*!

This is expected to be the residual equation for a Leaf PK, therefore it lives
on a domain.

It is expected to be of the form:


r = b_1 + ... b_k + rhs_A00 + rhs_A01 + ... rhs_Anm - (A00 x0 + A01 x0 + ...
A0n x0 + A11 x1 + ... Anm xn)

Where:

  x_0                  | is the primary variable associated with this residual,
                       | e.g. the diagonal entry in a coupled operator
  x_j, j=1...NJ        | are off-diagonal primary variables.
  A0i, rhs_A0i,        | are diagonal, local operators and their RHS boundary
    i=0...NI           | conditions/sources, i = 0...NI
  Ajm, rhs_Ajm         | are not-necessarily diagonal, local operators and their
    m=0...NM           | RHSs applied to the x_j, e.g. j = 1...NJ
  b_k, k=0...NK        | are arbitrary vector contributions (e.g. lumped mass
                       | matrix accumulation terms, source terms, which CANNOT
                       | NOT BE AFFECTED BY boundary conditions (i.e, for
                       | conservation equations discretized using control volume
                       | methods, BCs affect only faces while sources are on cells.

Notes:

- All of NJ, NI, NM, and NK are arbitrarily up to the user.

- The domain and range of the A0i must be subsets of the r space.  We take r's
  space to be the union of the domain and range spaces of the A0i and the
  spaces of the rhs_k.  We prefer to take the r space this way as it allows
  discretizations to impose the space -- for instance, in a conservation
  equation, FV methods might only need cells, while MFD methods need cells and
  faces.  This need not be known by the PK.

- The ranges of the Ajm (likewise the rhs_Ajm spaces) must be subsets of the r
space.

- The spaces of the b_k must be subsets of the r space.


Example 1:
~~~~~~~~~~
Standard parabolic PDE e.g. Richards equation.

Assume Richards equation is of the form:

.. math::
    \frac{\partial \Psi }{\partial t} - \nabla \cdot k(p) K \nabla p = Q

Where :math:`\Psi = \phi \rho s(p)` is the water content (the conserved
quantity), p is the pressure, k and K are the relative and absolutle
permeabilities, and Q is a volumetric source of water.

In this case:

- x0 is p, the pressure field,
- A00 is the only operator, the diffusion operator  :math:`\nabla \cdot k(p) K \nabla`,
- rhs_A00 accepts any flux boundary conditions for the operator,
- b_0 is the accumulation vector, dPsi/dt
- b_1 is the source Q


Example 2:
~~~~~~~~~~
Energy equation with an advection term driven by the Darcy flux.

Assume that energy transport is being coupled to Richards equation from Example
1.  In that case, we assume there are two primary variables -- temperature T is
the "diagonal" primary variable while pressure p is the off-diagonal primary
variable.  Then the equation is given by:

.. math::
    \frac{\partial E(p,T) }{\partial t} - \nabla \cdot \kappa \nabla T + \nabla \cdot h(T) q_{darcy} = Q_E + h(T_{source}) Q

In this case:

- x0 is T, the temperature field
- x1 is p, the pressure field
- A00 is the diagonal diffusion term, -\nabla \cdot \kappa \nabla T
- rhs_A00 includes Neumann energy diffusive flux boundary conditions
- A10 is the advective term, -\nabla \cdot h(T) k(p) K \nabla
- rhs_A10 includes advective flux boundary conditions (typically one rhs_A00 or rhs_A10 is 0)
- b_0 is the accumulation vector, dE/dt
- b_1 is the energy source, Q_E
- b_2 is the enthalpy associated with the mass source, h(T)Q

By splitting out terms in this formal sense, we can ask for derivatives of this
PDE.

- with respect to T: sum_k( db_k/dT ) + A00
- with respect to p: sum_k( db_k/dp ) + A10

*/


#ifndef STATE_EVALUATOR_OPERATOR_APPLY_HH_
#define STATE_EVALUATOR_OPERATOR_APPLY_HH_

#include "Debugger.hh"
#include "EvaluatorSecondary.hh"
#include "Evaluator_Factory.hh"

namespace Amanzi {

class Evaluator_OperatorApply : public EvaluatorSecondary {
 public:
  explicit Evaluator_OperatorApply(Teuchos::ParameterList& plist);
  Evaluator_OperatorApply(const Evaluator_OperatorApply& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new Evaluator_OperatorApply(*this));
  }

  virtual void EnsureCompatibility(State& S) override;

  // handled in ensure compatibility?
  //  virtual void EnsureCompatibleDerivative(State &S, const Key& wrt_key,
  //  const Key& wrt_tag) override {}

  virtual std::string name() const override { return "operator application"; }

  virtual bool UpdateDerivative(State& S, const Key& requestor,
          const Key& wrt_key, const Key& wrt_tag) override;

 protected:
  // These do the actual work
  virtual void Update_(State& S) override;

  // creates the operator for applying inverses
  virtual void
  UpdateDerivative_(State& S, const Key& wrt_key, const Key& wrt_tag) override;

 protected:
  Key global_operator_key_;          // key of operator for calculating the A_0i * x_0
  Key x0_key_;                       // x_0
  Teuchos::Array<Key> x_keys_;       // x_i
  Teuchos::Array<Key> op0_keys_;     // A_0i
  Teuchos::Array<Key> op0_rhs_keys_; // rhs_A_0i

  std::vector<Teuchos::Array<Key>> op_keys_;     // A_ji
  std::vector<Teuchos::Array<Key>> op_rhs_keys_; // rhs_A_ji
  Teuchos::Array<double> rhs_scalars_; // scale the rhs_A_ji (allows +1 for
                                       // sources, -1 for accumulation terms)

  Teuchos::Array<Key> rhs_keys_; // rhs_k

  std::string primary_entity_;
  AmanziMesh::Entity_kind primary_entity_kind_;

  // debugger for dumping vectors
  Teuchos::RCP<Debugger> db_;

private:
  static Utils::RegisteredFactory<Evaluator, Evaluator_OperatorApply> fac_;
};

} // namespace Amanzi

#endif
