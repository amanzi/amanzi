/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_PRECONDITIONER_BLOCK_ILU_HH_
#define AMANZI_PRECONDITIONER_BLOCK_ILU_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Factory.hpp"

#include "exceptions.hh"
#include "Preconditioner.hh"

#include "cuda_decl.h"

namespace Amanzi {
namespace AmanziSolvers {

using Ifpack2_PC_type = Ifpack2::Preconditioner<Matrix_type::scalar_type,
                                               Matrix_type::local_ordinal_type,
                                               Matrix_type::global_ordinal_type,
                                               Matrix_type::node_type>;

class PreconditionerIfpack2 : public Preconditioner {
 public:
  PreconditionerIfpack2():
    Preconditioner(), initialized_(false) {};

  virtual void set_inverse_parameters(Teuchos::ParameterList& plist) override {
    plist_ = plist;
    std::string vo_name = this->name()+" ("+plist_.get<std::string>("method")+")";
    vo_ = Teuchos::rcp(new VerboseObject(vo_name, plist_));
    initialized_ = true;
  }

  virtual void initializeInverse() override {
    nvtxRangePush("init-inv");
    Ifpack2::Factory factory;
    A_ = h_;
    std::string method = plist_.get<std::string>("method");
    pc_ = factory.create(method, A_);
    pc_->setParameters(plist_);
    pc_->initialize();
    if (vo_->os_OK(Teuchos::VERB_HIGH)) plist_.print(*vo_->os());
    nvtxRangePop(); 
  }

  virtual void computeInverse() override {
    nvtxRangePush("comp-inv"); 
    pc_->compute();
    if (vo_->os_OK(Teuchos::VERB_HIGH)) pc_->describe(*vo_->os(), vo_->getVerbLevel());
    nvtxRangePop(); 
  }

  virtual int returned_code() const override final { return returned_code_; }

  virtual int applyInverse(const Vector_type& v, Vector_type& hv) const override {
    nvtxRangePush("app-inv"); 

    pc_->apply(v, hv);
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) pc_->describe(*vo_->os(), vo_->getVerbLevel());
    nvtxRangePop(); 
    return 0;
  }

  virtual std::string returned_code_string() const override
  {
    switch (returned_code()) {
      case -1 :
        return "Generic Ifpack2 error.";
      case -2 :
        return "Ifpack2 says input data not valid.";
      case -3 :
        return "Ifpack2 says data not correctly preprocessed.";
      case -4 :
        return "Ifpack2 says problem encountered during algorithm, e.g. divide-by-zero, out-of-bounds, etc.";
      case -5 :
        return "Ifpack2 out-of-memory";
      case 0:
        return "Ifpack2 not yet applied.";
      case 1:
        return "success";
    }
    return "unknown error";
  }

 protected:
  Teuchos::ParameterList plist_;

  Teuchos::RCP<Ifpack2_PC_type> pc_;
  Teuchos::RCP<VerboseObject> vo_;

  Teuchos::RCP<const Matrix_type> A_;


  bool initialized_;
  mutable int returned_code_;

};

} // namespace AmanziPreconditioners
} // namespace Amanzi

#endif
