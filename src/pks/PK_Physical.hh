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

  PK_Physical(const Comm_ptr_type& comm,
              Teuchos::ParameterList& pk_tree,
              const Teuchos::RCP<Teuchos::ParameterList>& glist,
              const Teuchos::RCP<State>& S)
    : PK(comm, pk_tree, glist, S)
  {};

  // Virtual destructor
  virtual ~PK_Physical() = default;

  // Default implementations of PK methods.
  // -- transfer operators -- pointer copies only
  virtual Teuchos::RCP<TreeVectorSpace>
  getSolutionSpace() const override;

  virtual void State_to_Solution(const Tag& tag, TreeVector& soln) override;
  virtual void Solution_to_State(const TreeVector& soln, const Tag& tag) override;

  // access
  Key domain() { return domain_; }
  Teuchos::RCP<Debugger> debugger() { return db_; }

 protected:
  virtual void ParseParameterList_() override {
    domain_ = plist_->get<std::string>("domain name", "domain");
    mesh_ = S_->GetMesh(domain_);
  }

 protected:
  // name of domain, associated mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Key domain_;

  // solution and evaluator
  std::string key_;

  // debugger for dumping vectors
  Teuchos::RCP<Debugger> db_;
};

} // namespace Amanzi

#endif
