/*
  Process Kernels

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov, Ethan Coon

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

class PK_BDF : virtual public PK,
	       public BDFFnBase<TreeVector> {
 public:
  PK_BDF()
    : PK(),
      BDFFnBase<TreeVector>() {};

  PK_BDF(Teuchos::ParameterList& pk_tree,
	 const Teuchos::RCP<Teuchos::ParameterList>& glist,
	 const Teuchos::RCP<State>& S,
	 const Teuchos::RCP<TreeVector>& soln)
    : PK(pk_tree, glist, S, soln),
      BDFFnBase<TreeVector>() {};

  // access to operators and PDEs in sub-PKs
  virtual Teuchos::RCP<Operators::Operator>
      my_operator(const Operators::OperatorType& type) { return Teuchos::null; }

  virtual Teuchos::RCP<Operators::PDE_HelperDiscretization>
      my_pde(const Operators::PDEType& type) { return Teuchos::null; }

};

}  // namespace Amanzi

#endif
