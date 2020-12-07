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

class PreconditionerMueLu : public Preconditioner {
 public:
  PreconditionerMueLu() :
      Preconditioner(),
      initialized_(false) {};

  virtual ~PreconditionerMueLu() {}

  virtual void set_matrices(const Teuchos::RCP<Matrix_type>& m,
			    const Teuchos::RCP<Matrix_type>& h) override final;

  virtual void set_inverse_parameters(Teuchos::ParameterList& list) override final {
    list_ = list; 
    std::string vo_name = this->name()+" ("+plist_.get<std::string>("method")+")";
    vo_ = Teuchos::rcp(new VerboseObject(vo_name, plist_));
    initialized_ = true;
  }

  virtual void initializeInverse() override final{
    A_ = h_; 
    pc_ = MueLu::CreateTpetraPreconditioner(A_,plist_); 

  }
  virtual void computeInverse() override final{
    // pc_->compute(); ? 
  }
  virtual int applyInverse(const Vector_type& v, Vector_type& hv) const override final{
    pc_->apply(v, hv);
  }

  virtual int returned_code() const override final { return returned_code_; }
  virtual std::string returned_code_string() const override final {
    if (returned_code_ == 0) return "success";
    return "PreconditionerMueLu: unknown error";
  }

 private:

  using toperator = Tpetra::Operator<double,int,int>;
  using moperator = MueLu::TpetraOperator<double,int,int>; 

  Teuchos::ParameterList list_;
  Teuchos::RCP<moperator> pc_; 
  Teuchos::RCP<topetator> A_; 

  mutable int returned_code_;
  bool initialized_;
};

}  // namespace AmanziSolvers
}  // namespace Amanzi



#endif
