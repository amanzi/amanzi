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
#include "PDE_HelperDiscretization.hh"
#include "EvaluatorPrimary.hh"

#include "PK_Physical.hh"

namespace Amanzi {

class PK_BDF : virtual public PK, public BDFFnBase<TreeVector> {
 public:
  PK_BDF() : PK(), BDFFnBase<TreeVector>(){};

  PK_BDF(const Comm_ptr_type& comm,
         Teuchos::ParameterList& pk_tree,
         const Teuchos::RCP<Teuchos::ParameterList>& glist,
         const Teuchos::RCP<State>& S)
    : PK(comm, pk_tree, glist, S),
      BDFFnBase<TreeVector>()
    {};

  // access to operators and PDEs in sub-PKs
  virtual Teuchos::RCP<Operators::Operator> my_operator(const Operators::Operator_kind& type)
  {
    return Teuchos::null;
  }

  virtual Teuchos::RCP<Operators::PDE_HelperDiscretization> my_pde(const Operators::PDE_kind& type)
  {
    return Teuchos::null;
  }
};

} // namespace Amanzi

#endif
