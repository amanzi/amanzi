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

r = rhs_1 + ... rhs_k + rhs_A00 + rhs_A01 + ... rhs_Anm - (A00 x0 + A01 x0 + ...
A0n x0 + A11 x1 + ... Anm xn)

Where:

  x0                   | is the primary variable associated with this residual
  x1 ... xn            | are off-diagonal primary variables.
  A0i, rhs_A0i         | are diagonal, local operators and their RHS boundary
                       | conditions/sources.
  Aji, rhs_Aji         | are not-necessarily diagonal, local operators and their
                       | RHSs applied to the xj.
  rhs_1 ... rhs_k      | are arbitrary source terms which CANNOT NOT BE AFFECTED
BY | boundary conditions (i.e, for conservation equations | discretized using
control volume methods, BCs affect only | faces while sources are on cells.

Note that we can infer some constraints here:

- The domain and range of the A0i must be subsets of the r space.  We take r's
  space to be the union of the domain and range spaces of the A0i and the
  spaces of the rhs_k.  We prefer to take the r space this way as it allows
  discretizations to impose the space -- for instance, in a conservation
  equation, FV methods might only need cells, while MFD methods need cells and
  faces.  This need not be known by the PK.

- The ranges of the Aji (likewise the rhs_Aji spaces) must be subsets of the r
space.

- The spaces of the rhs_k must be subsets of the r space.



NOTE: It is unclear what to do about transformed primary variables at this
point, i.e. overland flow which writes operator(h(p)), and therefore needs to
scale the linear operators by dh/dp. --etc

*/


#ifndef STATE_EVALUATOR_OPERATOR_APPLY_HH_
#define STATE_EVALUATOR_OPERATOR_APPLY_HH_

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

 protected:
  // These do the actual work
  virtual void Update_(State& S) override;

  // creates the operator for applying inverses
  virtual void
  UpdateDerivative_(State& S, const Key& wrt_key, const Key& wrt_tag) override;

 protected:
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

 private:
  static Utils::RegisteredFactory<Evaluator, Evaluator_OperatorApply> fac_;
};

} // namespace Amanzi

#endif
