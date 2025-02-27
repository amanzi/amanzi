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

namespace Amanzi {

class PK_Physical : virtual public PK {
 public:
  PK_Physical() : PK(){};

  PK_Physical(Teuchos::ParameterList& pk_tree,
              const Teuchos::RCP<Teuchos::ParameterList>& glist,
              const Teuchos::RCP<State>& S,
              const Teuchos::RCP<TreeVector>& soln)
    : PK(pk_tree, glist, S, soln)
  {
    domain_ = Keys::readDomain(*plist_, "domain", "domain");
    mesh_ = S_->GetMesh(domain_);
  };

  // Virtual destructor
  virtual ~PK_Physical() = default;

  // Default implementations of PK methods.
  // -- transfer operators -- pointer copies only
  virtual void State_to_Solution(const Tag& tag, TreeVector& soln) override;
  virtual void Solution_to_State(const TreeVector& soln, const Tag& tag) override;

  virtual void parseParameterList() override;
  virtual void Setup() override;
  virtual void Initialize() override;


  virtual void CommitStep(double t_old, double t_new, const Tag& tag_next) override;
  virtual void FailStep(double t_old, double t_new, const Tag& tag_next) override;

  // access
  Key domain() { return domain_; }
  Teuchos::RCP<Debugger> debugger() { return db_; }

 protected:
  void ChangedSolutionPK_(const Tag& tag);

  // name of domain, associated mesh
  Key domain_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  // solution and evaluator
  Key key_;
  Key passwd_;

  // debugger for dumping vectors
  Teuchos::RCP<Debugger> db_;
};


} // namespace Amanzi

#endif
