/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_PRECONDITIONER_IDENTITY_HH_
#define AMANZI_PRECONDITIONER_IDENTITY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

template <class Matrix, class Vector>
class PreconditionerIdentity : public Preconditioner<Matrix, Vector> {
 public:
  PreconditionerIdentity(){};
  ~PreconditionerIdentity(){};

  void
  Init(const std::string& name, const Teuchos::ParameterList& list) override{};
  void Update(const Teuchos::RCP<Matrix>& A) override{};
  void Destroy() override{};

  int applyInverse(const Vector& v, Vector& hv) const override
  {
    hv.assign(v);
    return 0;
  }

  int returned_code() override { return 0; }
};


} // namespace AmanziPreconditioners
} // namespace Amanzi

#endif
