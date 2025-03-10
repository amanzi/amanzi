/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Process Kernels

  Default base with a few methods implemented in standard ways.
*/

#ifndef AMANZI_PK_PHYSICAL_HH_
#define AMANZI_PK_PHYSICAL_HH_

#include <string>

#include "Teuchos_ParameterList.hpp"
#include "VerboseObject.hh"

#include "Debugger.hh"
#include "Key.hh"
#include "EvaluatorPrimary.hh"
#include "PK.hh"
#include "PK_Physical.hh"

namespace Amanzi {

class PK_Physical_Default : public PK_Physical {
 public:
  PK_Physical_Default(Teuchos::ParameterList& pk_tree,
                      const Teuchos::RCP<Teuchos::ParameterList>& glist,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& soln)
    : PK_Physical(pk_tree, glist, S, soln),
      PK(pk_tree, glist, S, soln)
  {};

  // Default implementations of PK methods.
  virtual void parseParameterList() override;
  virtual void Setup() override;
  virtual void Initialize() override;

  virtual void CommitStep(double t_old, double t_new, const Tag& tag_next) override;
  virtual void FailStep(double t_old, double t_new, const Tag& tag_next) override;

  void ChangedSolutionPK(const Tag& tag) override;

 protected:
  // solution and evaluator
  Key passwd_;
};


} // namespace Amanzi

#endif
