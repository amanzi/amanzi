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

#pragma once

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

class PK_Default : public PK {
 public:
  PK_Default() {};
  // Required constructor for use by the PK factory.
  PK_Default(const Comm_ptr_type& comm,
             Teuchos::ParameterList& pk_tree,
             const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
             const Teuchos::RCP<State>& S);

  // Virtual destructor
  virtual ~PK() = default;

  // Return PK's name
  virtual const std::string& getName() const { return name_; }
  virtual const std::string& getType() const = 0;

  // Set a tag interval for advancing
  // Set the tags to integrate between
  virtual void setTags(const Tag& current, const Tag& next)
  {
    tag_current_ = current;
    tag_next_ = next;
  }

 protected:
  Comm_ptr_type comm_;
  Teuchos::RCP<Teuchos::ParameterList> plist_;
  std::string name_;
  Tag tag_current_, tag_next_; // tags for time integration

  Teuchos::RCP<State> S_;          // global data manager
  Teuchos::RCP<VerboseObject> vo_; // fancy IO
};

} // namespace Amanzi

