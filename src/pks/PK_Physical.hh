/*
  Process Kernels

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

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
  PK_Physical()
    : PK() {};

  PK_Physical(Teuchos::ParameterList& pk_tree,
              const Teuchos::RCP<Teuchos::ParameterList>& glist,
              const Teuchos::RCP<State>& S,
              const Teuchos::RCP<TreeVector>& soln)
    : PK(pk_tree, glist, S, soln)
  {
    domain_ = plist_->get<std::string>("domain name", "domain");
    mesh_ = S_->GetMesh(domain_);
  };

  // Virtual destructor
  virtual ~PK_Physical() = default;

  // Default implementations of PK methods.
  // -- transfer operators -- pointer copies only
  virtual void State_to_Solution(const Tag& tag, TreeVector& soln) override;
  virtual void Solution_to_State(const TreeVector& soln, const Tag& tag) override;

  // access
  Key domain() { return domain_; }
  Teuchos::RCP<Debugger> debugger() { return db_; }

 protected:
  // Helper method to add a primary variable evaluator
  void AddDefaultPrimaryEvaluator_(const Key& key, const Tag& tag = Tags::DEFAULT);
  void AddDefaultPrimaryEvaluator_(const Tag& tag = Tags::DEFAULT) { AddDefaultPrimaryEvaluator_(key_, tag); }

protected:
  // name of domain, associated mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Key domain_;

  // solution and evaluator
  std::string key_;

  // debugger for dumping vectors
  Teuchos::RCP<Debugger> db_;
};


// non-meber Helper method to initialize a CV field
void InitializeCVField(const Teuchos::RCP<State>& S, const VerboseObject& vo,
                       const Key& key, const Tag& tag, const Key& passwd,
                       double default_val);

} // namespace Amanzi

#endif
