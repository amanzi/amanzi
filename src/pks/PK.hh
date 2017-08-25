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

``PK`` is a virtual interface for a Process Kernel. Note that PKs
  deriving from this class must implement the constructor
  interface as well, and should add the private static member
  (following the Usage notes in src/pks/PK_Factory.hh) to register the
  derived PK with the PK factory.

*/

#ifndef AMANZI_PK_HH_
#define AMANZI_PK_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Key.hh"
#include "BDFFnBase.hh"
#include "Explicit_TI_FnBase.hh"

namespace Amanzi {

class Debugger;
class TreeVector;

class PK {
 public:

  // Virtual destructor
  virtual ~PK() = default;

  // construct all sub-PKs.  This is not a part of the constructor as it must
  // be virtual.
  virtual void ConstructChildren() = 0;
  
  // Setup: forms the DAG, pushes meta-data into State
  virtual void Setup() = 0;
  
  // Initialize: set values for owned variables.
  virtual void Initialize() = 0;

  // Advance PK from time tag old to time tag new
  virtual bool AdvanceStep(const Key& tag_old, const Key& tag_new) = 0;

  // Returns validity of the step taken from tag_old to tag_new
  virtual bool ValidStep(const Key& tag_old, const Key& tag_new) = 0;

  // Do work that can only be done if we know the step was successful.
  virtual void CommitStep(const Key& tag_old, const Key& tag_new) = 0;

  // Revert a step from tag_new back to tag_old
  virtual void FailStep(const Key& tag_old, const Key& tag_new) = 0;

  // Choose a time step compatible with physics.
  virtual double get_dt() = 0;
  
  // Calculate any diagnostics at tag, currently used for visualization.
  virtual void CalculateDiagnostics(const Key& tag) = 0;

  // Mark, as changed, any primary variable evaluator owned by this PK
  virtual void ChangedSolutionPK(const Key& tag) = 0;

  virtual void StateToSolution(TreeVector& soln, const Key& tag) = 0;
  virtual void SolutionToState(TreeVector& soln, const Key& tag) = 0;
  
  // Return PK's name
  virtual std::string name() = 0;

  // Accessor for debugger, for use by coupling MPCs
  virtual Teuchos::Ptr<Debugger> debugger() = 0;

  
};

template<typename Vector=TreeVector>
class PK_BDF : public PK,
               public BDFFnBase<Vector> {
 public:
  virtual ~PK_BDF() = default;
};

template<typename Vector=TreeVector>
class PK_Explicit : public PK,
                    public Explicit_TI::fnBase<Vector> {
 public:
  virtual ~PK_Explicit() = default;
};
  




}  // namespace Amanzi

#endif
