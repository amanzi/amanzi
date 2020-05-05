/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! A preconditioner that assembles, then calls an assembled preconditioner.

#ifndef AMANZI_PRECONDITIONER_ASSEMBLED_HH_
#define AMANZI_PRECONDITIONER_ASSEMBLED_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Tpetra_RowMatrix.hpp"

#include "AmanziTypes.hh"
#include "SuperMap.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

template <class Matrix, class Vector>
class PreconditionerAssembled : public Preconditioner<Matrix, Vector> {
 public:

  PreconditionerAssembled() {};
  ~PreconditionerAssembled() {};

  void
  Init(const std::string& name,
       const Teuchos::RCP<Teuchos::ParameterList>& plist) override {
    PreconditionerFactory<Matrix_type,Vector_type> fac;
    pc_ = fac.Create(plist);
  };

  void Update(const Teuchos::RCP<Matrix>& mat) override
  {
    mat->AssembleMatrix();
    A_ = mat->A();
    pc_->Update(A_);
    smap_ = mat->getSuperMap();
  }

  void Destroy() override {};
  int applyInverse(const Vector& Y, Vector& X) const override
  {
    Vector_type Xcopy(A_->getRowMap());
    Vector_type Ycopy(A_->getRowMap());

    returned_code_ = copyToSuperVector(*smap_, Y, Ycopy);
    returned_code_ |= pc_->applyInverse(Ycopy, Xcopy);
    returned_code_ |= copyFromSuperVector(*smap_, Xcopy, X);
    return returned_code_;
  }

  int returned_code() override { return returned_code_; }

 private:
  Teuchos::RCP<Preconditioner<Matrix_type,Vector_type>> pc_;
  Teuchos::RCP<Matrix_type> A_;
  Teuchos::RCP<const Operators::SuperMap> smap_;
  mutable int returned_code_;
};


template<class Matrix>
class PreconditionerAssembled<Matrix, Vector_type>
    : public Preconditioner<Matrix, Vector_type> {
 public:

  PreconditionerAssembled() {};
  ~PreconditionerAssembled() {};

  void
  Init(const std::string& name,
       const Teuchos::RCP<Teuchos::ParameterList>& plist) override {
    PreconditionerFactory<Matrix_type,Vector_type> fac;
    pc_ = fac.Create(plist);
  };

  void Update(const Teuchos::RCP<Matrix>& mat) override
  {
    mat->AssembleMatrix();
    pc_->Update(mat->A());
  }

  void Destroy() override {};

  int applyInverse(const Vector_type& Y, Vector_type& X) const override
  {
    returned_code_ |= pc_->applyInverse(Y, X);
    return returned_code_;
  }

  int returned_code() override { return returned_code_; }

 private:
  Teuchos::RCP<Preconditioner<Matrix_type,Vector_type>> pc_;
  mutable int returned_code_;
};
  


} // namespace AmanziPreconditioners
} // namespace Amanzi

#endif
