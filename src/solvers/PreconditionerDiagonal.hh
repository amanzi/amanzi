/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_PRECONDITIONER_DIAGONAL_HH_
#define AMANZI_PRECONDITIONER_DIAGONAL_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

template <class Matrix, class Vector>
class PreconditionerDiagonal : public Preconditioner<Matrix, Vector> {
 public:
  PreconditionerDiagonal(){};
  ~PreconditionerDiagonal(){};

  void
  Init(const std::string& name, const Teuchos::RCP<Teuchos::ParameterList>& plist) override{};
  void Update(const Teuchos::RCP<Matrix>& A) override
  {
    work_vec_ = Teuchos::rcp(new Vector(A->getRowMap()));
    work_vec_->putScalar(0.);
    A->getLocalDiagCopy(*work_vec_);
    work_vec_->reciprocal(*work_vec_);
  }

  void Destroy() override{};
  int applyInverse(const Vector& v, Vector& hv) const override
  {
    hv.elementWiseMultiply(1., v, *work_vec_, 0.);
    return 0;
  }
  int returned_code() override { return 0; }

 private:
  Teuchos::RCP<Vector> work_vec_;
};

} // namespace AmanziPreconditioners
} // namespace Amanzi

#endif
