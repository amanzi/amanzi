/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! The interface for a Process Kernel, an equation or system of equations.
/*!

A process kernel represents a single or system of partial/ordinary
differential equation(s) or conservation law(s), and is used as the
fundamental unit for coupling strategies.

Implementations of this interface typically are either an MPC
(multi-process coupler) whose job is to heirarchically couple several
other PKs and represent the system of equations, or a Physical PK,
which represents a single equation.

Note there are two PK specs -- the first is the "typed" spec, which appears in
the "cycle driver" list in the PK tree and has no other parameters other than
its type and its children.  The second is the spec for the base class PK, which
is inherited and included by each actual PK, lives in the "PKs" sublist of
"main", and has all needed parameters.

.. _pk-typed-spec:
.. admonition:: pk-typed-spec

    * `"PK type`" ``[string]`` One of the registered PK types

Example:

.. code-block:: xml

  <ParameterList name="PK tree">
    <ParameterList name="Top level MPC">
      <Parameter name="PK type" type="string" value="strong MPC"/>
      <ParameterList name="sub PK 1">
        ...
      </ParameterList>
      <ParameterList name="sub PK 2">
        ...
      </ParameterList>
      ...
    </ParameterList>
  </ParameterList>


.. _pk-spec:
.. admonition:: pk-spec

    * `"PK type`" ``[string]`` One of the registered PK types.  Note this must
      match the corresponding entry in the ``[pk-typed-spec]``
    * `"verbose object`" ``[verbose-object-spec]`` **optional** See `Verbose Object`_

Example:

.. code-block:: xml

  <ParameterList name="PKs">
    <ParameterList name="my cool PK">
      <Parameter name="PK type" type="string" value="my cool PK"/>
       ...
    </ParameterList>
  </ParameterList>

*/


/*
Developer's note:

``PK`` is a virtual interface for a Process Kernel. Note that PKs
  deriving from this class must implement the commented constructor
  interface as well, and should add the private static member
  (following the Usage notes in src/pks/PK_Factory.hh) to register the
  derived PK with the PK factory.
*/

#ifndef AMANZI_PK_HH_
#define AMANZI_PK_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Key.hh"
#include "Tag.hh"
#include "TreeVectorSpace.hh"
#include "StateDefs.hh"

namespace Amanzi {

class State;
class PK_Factory;
class TreeVector;

class PK {
 public:
  PK(){};
  // Required constructor for use by the PK factory.
  PK(const Comm_ptr_type& comm,
     Teuchos::ParameterList& pk_tree,
     const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
     const Teuchos::RCP<State>& S)
    : comm_(comm),
      name_(Keys::cleanPListName(pk_tree)),
      tag_current_(Tags::CURRENT),
      tag_next_(Tags::NEXT),
      S_(S)
  {
    Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(global_plist, "PKs", true);
    if (pk_list->isSublist(name_)) {
      plist_ = Teuchos::sublist(pk_list, name_);
    } else {
      Errors::Message msg;
      msg << "There is no sublist for PK \"" << name_ << "\" in PKs list";
      Exceptions::amanzi_throw(msg);
    }

    // set up the VerboseObject
    // Note, this allows for overriding the vo plist for individual PKs in a
    // collection of PKs
    Teuchos::RCP<Teuchos::ParameterList> vo_plist = plist_;
    if (plist_->isSublist(name_ + " verbose object")) {
      vo_plist = Teuchos::rcp(new Teuchos::ParameterList(*plist_));
      vo_plist->set("verbose object", plist_->sublist(name_ + " verbose object"));
    }

    //  some tests provide nullptr
    vo_ = Teuchos::rcp(new VerboseObject(comm_, name_, *vo_plist));
  };

  // Virtual destructor
  virtual ~PK() = default;

  // Setup, requiring all data used in state
  virtual void Setup() = 0;

  // Initialize owned (dependent) variables.
  virtual void Initialize() = 0;

  // Return PK's name
  virtual const std::string& getName() const { return name_; }
  virtual const std::string& getType() const = 0;

  // Choose a time step compatible with physics.
  virtual double getDt() = 0;

  // Set a time step for a PK.
  virtual void setDt(double dt) = 0;

  // Set a tag interval for advancing
  // Set the tags to integrate between
  virtual void setTags(const Tag& current, const Tag& next)
  {
    tag_current_ = current;
    tag_next_ = next;
  }

  // Advance PK from time t_old to time t_new. True value of the last
  // parameter indicates drastic change of boundary and/or source terms
  // that may need PK's attention.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) = 0;

  // This is called after ALL PKs have successfully advanced their
  // steps, so information needed to back up can be overwritten.
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) = 0;

  // This is called if ANY PK has failed; do what is needed to back up for a
  // new attempt at the step.
  virtual void FailStep(double t_old, double t_new, const Tag& tag) {}

  // Calculate any diagnostics at S->time(), currently for visualization.
  virtual void CalculateDiagnostics(const Tag& tag) {}

  /////////////////////////////////////////////////////////////////////
  // -- transfer operators
  virtual Teuchos::RCP<TreeVectorSpace>
  getSolutionSpace() const = 0;
  virtual void State_to_Solution(const Tag& tag, TreeVector& soln) = 0;
  virtual void Solution_to_State(const TreeVector& soln, const Tag& tag) = 0;

  // Tag the primary variable as changed in the DAG
  virtual void ChangedSolutionPK(const Tag& tag) { AMANZI_ASSERT(false); }

  // When including ValidStep() in Advance(), make this protected!  refs
  // amanzi/ats#110
  // Check whether the solution calculated for the new step is valid.
  virtual bool ValidStep() { return true; }

  virtual void ParseParameterList_() {}

 protected:
  Comm_ptr_type comm_;
  Teuchos::RCP<Teuchos::ParameterList> plist_;
  std::string name_;
  Tag tag_current_, tag_next_; // tags for time integration

  Teuchos::RCP<State> S_;             // global data manager
  Teuchos::RCP<VerboseObject> vo_; // fancy IO

};

} // namespace Amanzi

#endif
