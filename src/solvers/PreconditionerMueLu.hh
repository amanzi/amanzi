/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef AMANZI_PRECONDITIONER_MUELU_HH_
#define AMANZI_PRECONDITIONER_MUELU_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "MueLu.hpp"
#include "MueLu_TpetraOperator.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziSolvers {

using TpetraOperator_type = Tpetra::Operator<Matrix_type::scalar_type,
                                             Matrix_type::global_ordinal_type,
                                             Matrix_type::local_ordinal_type>;
using MueLuOperator_type = MueLu::TpetraOperator<Matrix_type::scalar_type,
                                                 Matrix_type::global_ordinal_type,
                                                 Matrix_type::local_ordinal_type>;

class PreconditionerMueLu : public Preconditioner {
 public:
  PreconditionerMueLu() : Preconditioner() {}

  virtual void set_inverse_parameters(Teuchos::ParameterList& plist) override final
  {
    plist_ = plist;
    std::string vo_name = this->name() + " (" + plist_.get<std::string>("method") + ")";
  }

  virtual void initializeInverse() override final
  {
    plist_.remove("verbose object");
    plist_.remove("method");
  }

  virtual void computeInverse() override final
  {
    Teuchos::RCP<TpetraOperator_type> t_op = h_;
    pc_ = MueLu::CreateTpetraPreconditioner(t_op, plist_);
  }

  virtual int applyInverse(const Vector_type& v, Vector_type& hv) const override final
  {
    pc_->apply(v, hv);
    returned_code_ = 0;
    return returned_code_;
  }

  virtual int returned_code() const override final { return returned_code_; }

  virtual std::string returned_code_string() const override final
  {
    if (returned_code_ == 0) return "success";
    return "PreconditionerMueLu: unknown error";
  }

 private:
  Teuchos::ParameterList plist_;
  //  Teuchos::RCP<MueLuOperator_type> pc_;
  Teuchos::RCP<TpetraOperator_type> pc_;


  mutable int returned_code_;
};

} // namespace AmanziSolvers
} // namespace Amanzi


#endif
