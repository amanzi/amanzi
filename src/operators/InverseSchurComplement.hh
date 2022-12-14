/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! Provides ApplyInverse() using a Schur complement.
#pragma once

#include "Inverse.hh"
#include "Operator.hh"

namespace Amanzi {
namespace AmanziSolvers {

//
// Class for assembled inverse methods.
//
class InverseSchurComplement : public Inverse<Operators::Operator,
                                              Operators::Operator,
                                              CompositeVector,
                                              CompositeVectorSpace> {
 public:
  InverseSchurComplement() {}

  virtual void set_inverse_parameters(Teuchos::ParameterList& plist) override final;

  virtual void InitializeInverse() override final;
  virtual void ComputeInverse() override final;
  virtual int ApplyInverse(const CompositeVector& y, CompositeVector& x) const override final;

  virtual double residual() const override final { return solver_->residual(); }

  virtual int num_itrs() const override final { return solver_->num_itrs(); }

  virtual void add_criteria(int criteria) override final { return solver_->add_criteria(criteria); }

  virtual int returned_code() const override final { return solver_->returned_code(); }

  virtual std::string returned_code_string() const override final
  {
    return solver_->returned_code_string();
  }

 protected:
  using Inverse<Operators::Operator, Operators::Operator, CompositeVector, CompositeVectorSpace>::
    h_;

  Teuchos::RCP<Inverse<Epetra_CrsMatrix, Epetra_CrsMatrix, Epetra_Vector, Epetra_Map>> solver_;
};


} // namespace AmanziSolvers
} // namespace Amanzi
