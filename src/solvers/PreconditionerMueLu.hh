/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#ifndef AMANZI_PRECONDITIONER_MUELU_HH_
#define AMANZI_PRECONDITIONER_MUELU_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "MueLu.hpp"
#include "MueLu_EpetraOperator.hpp"
#include "MueLu_CreateEpetraPreconditioner.hpp"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziSolvers {

class PreconditionerMueLu : public Preconditioner {
 public:
  PreconditionerMueLu() : Preconditioner() {};

  virtual void set_matrices(const Teuchos::RCP<Epetra_CrsMatrix>& m,
			    const Teuchos::RCP<Epetra_CrsMatrix>& h) override final;

  virtual void set_inverse_parameters(Teuchos::ParameterList& plist) override final {
    plist_ = plist;
    // std::string vo_name = this->name()+" ("+plist_.get<std::string>("method")+")";
  }

  virtual void InitializeInverse() override final;
  virtual void ComputeInverse() override final;
  virtual int ApplyInverse(const Vector_t& v, Vector_t& hv) const override final;

  virtual int returned_code() const override final { return returned_code_; }

  virtual std::string returned_code_string() const override final {
    if (returned_code_ == 0) return "success";
    return "PreconditionerMueLu: unknown error";
  }

 private:
  Teuchos::ParameterList plist_;
#if defined(HAVE_MUELU_EPETRA)
  Teuchos::RCP<MueLu::EpetraOperator> MueLu_;
#endif

  mutable int returned_code_;
};

}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
