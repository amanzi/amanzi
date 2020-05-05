/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Smoothers and incomplete factorizations from Ifpack2

#pragma once

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Factory.hpp"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

using Ifpack_PC_type = Ifpack2::Preconditioner<Matrix_type::scalar_type,
                                               Matrix_type::local_ordinal_type,
                                               Matrix_type::global_ordinal_type,
                                               Matrix_type::node_type>;

class PreconditionerIfpack2 : public Preconditioner<Matrix_type, Vector_type> {
 public:
  PreconditionerIfpack2() {};

  void
  Init(const std::string& name, const ParameterList_ptr_type& plist) override {
    plist_ = plist;

    if (Keys::startsWith(name, "ifpack2: ")) {
      name_ = name.substr(9,name.size());
    } else {
      name_ = name;
    }
  }

  void Update(const Teuchos::RCP<Matrix_type>& A) override {
    A_ = A;
    Ifpack2::Factory factory;
    pc_ = factory.create(name_, A_);
    pc_->setParameters(plist_->sublist(std::string("ifpack2: ")+name_+" parameters"));
    pc_->initialize();
    pc_->compute();
  }

  void Destroy() override {};

  int applyInverse(const Vector_type& v, Vector_type& hv) const override {
    pc_->apply(v, hv);
    return 0;
  }
    
  int returned_code() override { return 0; }

 protected:
  ParameterList_ptr_type plist_;
  std::string name_;

  Teuchos::RCP<Ifpack_PC_type> pc_;
  Teuchos::RCP<const Matrix_type> A_;
};

} // namespace AmanziPreconditioners
} // namespace Amanzi
