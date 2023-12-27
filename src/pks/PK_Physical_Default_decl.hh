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

#pragma once

#include <string>

#include "Teuchos_ParameterList.hpp"
#include "VerboseObject.hh"

#include "Key.hh"
#include "PK_Default.hh"

namespace Amanzi {

// forward declarations
class Debugger;
namespace Operators {
class BCs;
class PDE_HelperDiscretization;
enum class PDE_kind;
} // namespace Operators

//
// Note, this could be:
//  - PK_Default
template<class PK_type>
class PK_Physical_Default : public PK_type {
 public:
  PK_Physical_Default() : PK_type() {};

  PK_Physical_Default(const Comm_ptr_type& comm,
              Teuchos::ParameterList& pk_tree,
              const Teuchos::RCP<Teuchos::ParameterList>& glist,
              const Teuchos::RCP<State>& S)
    : PK_type(comm, pk_tree, glist, S) {};

  // Default implementations of PK methods.
  // -- transfer operators -- pointer copies only
  virtual Teuchos::RCP<TreeVectorSpace> getSolutionSpace() const override;
  virtual void moveStateToSolution(const Tag& tag, TreeVector& soln) override;
  virtual void moveSolutionToState(const TreeVector& soln, const Tag& tag) override;

  // access
  Key getDomain() { return domain_; }
  Teuchos::RCP<Debugger> getDebugger() { return db_; }
  Teuchos::RCP<Operators::BCs> getBCs() { return bc_; }

  virtual Teuchos::RCP<Operators::PDE_HelperDiscretization> getPDE(const Operators::PDE_kind& type) {
    return Teuchos::null;
  }

  virtual void parseParameterList() override;

  virtual bool isValidStep() override;

  // Tag the primary variable as changed in the DAG
  virtual void markChangedSolutionPK(const Tag& tag) override;

  virtual void setup() override;
  virtual void initialize() override;
  virtual void commitStep(double t_old, double t_new, const Tag& tag) override;
  virtual void failStep(double t_old, double t_new, const Tag& tag) override;

 protected:
  using PK_type::S_;
  using PK_type::tag_next_;
  using PK_type::tag_current_;
  using PK_type::name_;
  using PK_type::vo_;
  using PK_type::plist_;

  // name of domain, associated mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Key domain_;

  // primary variable
  std::string key_;

  // debugger for dumping vectors
  Teuchos::RCP<Debugger> db_;

  // step validity -- limit the max change of key_ in a given step
  double max_valid_change_;

  // boundary conditions
  Teuchos::RCP<Operators::BCs> bc_;
};

} // namespace Amanzi

