/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_EPETRAMATRIXDEFAULT_HH_
#define AMANZI_EPETRAMATRIXDEFAULT_HH_

#include "Teuchos_RCP.hpp"
#include "errors.hh"
#include "EpetraMatrix.hh"

namespace Amanzi {

template <class Matrix>
class EpetraMatrixDefault : public EpetraMatrix {
 public:
  // constructor
  EpetraMatrixDefault(const Teuchos::ParameterList& plist) : plist_(plist)
  {
    AmanziPreconditioners::PreconditionerFactory fac;
    if (!plist.isSublist("preconditioner")) {
      Errors::Message msg(
        "EpetraMatrixDefault: missing \"preconditioner\" sublist");
      Exceptions::amanzi_throw(msg);
    }
    pc_ = fac.Create(plist_.sublist("preconditioner"));
  }

  // copy constructor
  EpetraMatrixDefault(const EpetraMatrixDefault<Matrix>& other)
    : plist_(other.plist_)
  {
    AmanziPreconditioners::PreconditionerFactory fac;
    pc_ = fac.Create(plist_);
    if (other.m_ != Teuchos::null) {
      Teuchos::RCP<Matrix> m_copy = Teuchos::rcp(new Matrix(*other.m_));
      Update(m_copy);
    }
  }

  // Methods to work with the PC
  void Update(const Teuchos::RCP<Matrix>& m)
  {
    m_ = m;
    pc_->Update(m);
  }

  void Destroy() { pc_->Destroy(); }

  // Interface to work with LinearOperator

  // Vector space of the Matrix's domain.
  virtual const Epetra_BlockMap& DomainMap() const { return m_->DomainMap(); }

  // Vector space of the Matrix's range.
  virtual const Epetra_BlockMap& RangeMap() const { return m_->RangeMap(); }

  // Virtual copy constructor.
  virtual Teuchos::RCP<EpetraMatrix> Clone() const
  {
    return Teuchos::rcp(new EpetraMatrixDefault<Matrix>(*this));
  }

  // Apply matrix, b <-- Ax, returns ierr
  virtual int Apply(const Epetra_Vector& x, Epetra_Vector& b) const
  {
    return m_->Apply(x, b);
  }

  // Apply the inverse, x <-- A^-1 b, returns ierr
  virtual int ApplyInverse(const Epetra_Vector& b, Epetra_Vector& x) const
  {
    return pc_->ApplyInverse(b, x);
  }

 protected:
  Teuchos::RCP<AmanziPreconditioners::Preconditioner> pc_;
  Teuchos::RCP<Matrix> m_;
  Teuchos::ParameterList plist_;
};

} // namespace Amanzi

#endif
