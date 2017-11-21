/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)
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

All PKs have the following parameters in their spec:

* `"PK type`" ``[string]``

  The PK type is a special key-word which corresponds to a given class in the PK factory.  See available PK types listed below.

* `"PK name`" ``[string]`` **LIST-NAME**

  This is automatically written as the `"name`" attribute of the containing PK sublist, and need not be included by the user.

Example:

.. code-block:: xml

  <ParameterList name="PKs">
    <ParameterList name="my cool PK">
      <Parameter name="PK type" type="string" value="my cool PK"/>
       ...
    </ParameterList>
  </ParameterList>

.. code-block:: xml

  <ParameterList name="PKs">
    <ParameterList name="Top level MPC">
      <Parameter name="PK type" type="string" value="strong MPC"/>
       ...
    </ParameterList>
  </ParameterList>

 */


/*
Developer's note:

This sets defaults for all things that can reasonably be ignored by other
mixins, and is typically the base-most class.
*/

#ifndef AMANZI_PK_DEFAULT_HH_
#define AMANZI_PK_DEFAULT_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "StateDefs.hh"
#include "TreeVector.hh"

namespace Amanzi {

class State;
class Debugger;

class PK_Default {
 public:

  // lone constructor
  PK_Default(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
     const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
     const Teuchos::RCP<State>& S,
     const Teuchos::RCP<TreeVector>& solution);

  // Setup: forms the DAG, pushes meta-data into State
  void Setup(const TreeVector& soln) {}
  
  // Initialize: set values for owned variables.
  void Initialize() {}
  
  // Returns validity of the step taken from tag_old to tag_new
  bool ValidStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                 const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new) { return true; }

  // Do work that can only be done if we know the step was successful.
  void CommitStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                  const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new) {}

  // Revert a step from tag_new back to tag_old
  void FailStep(const Key& tag_old, const Teuchos::RCP<TreeVector>& soln_old,
                const Key& tag_new, const Teuchos::RCP<TreeVector>& soln_new) {}

  // Calculate any diagnostics at tag, currently used for visualization.
  void CalculateDiagnostics(const Key& tag) {}

  // Return PK's name
  std::string name() { return name_; }

  // Accessor for debugger, for use by coupling MPCs
  Teuchos::Ptr<Debugger> debugger() { return Teuchos::null; }
  
 protected:
  // my subtree of the solution vector
  //  Teuchos::RCP<TreeVector> solution_;

  // state
  Teuchos::RCP<State> S_;

  // fancy IO
  Teuchos::RCP<VerboseObject> vo_;

  // my parameterlist
  Teuchos::RCP<Teuchos::ParameterList> plist_;
  
  // my name
  std::string name_;
  
};

}  // namespace Amanzi

#endif
