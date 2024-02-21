/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov, Ethan Coon
*/

/*
  Process Kernels

  This is a purely virtual base class for process kernels which use
  (implicit) time integrators.
*/

#ifndef AMANZI_PK_BDF_HH_
#define AMANZI_PK_BDF_HH_

#include "Teuchos_RCP.hpp"

#include "BDFFnBase.hh"
#include "Key.hh"
#include "Operator.hh"
#include "OperatorDefs.hh"

namespace Amanzi {

class PK_BDF : public PK, public BDFFnBase<TreeVector> {
 public:

  // access to operators and PDEs in sub-PKs
  virtual Teuchos::RCP<Operators::Operator> getOperator(const Operators::Operator_kind& type) = 0;
};

} // namespace Amanzi

#endif
