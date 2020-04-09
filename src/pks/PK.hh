/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! The interface for a Process Kernel, an equation or system of equations.

/*!
A process kernel represents a single or system of partial/ordinary
differential equation(s) or conservation law(s), and is used as the
fundamental unit for coupling strategies.

Three purely-virtual interfaces are defined here:

 - ``PK``, the default purely-virtual interface
 - ``PK_Explicit``, the composition of PK and Explicit_TI::fnBase
 - ``PK_Implicit``, the composition of PK and BDFFnBase


Implementations of this interface typically are either an MPC
(multi-process coupler) whose job is to heirarchically couple several
other PKs and represent the system of equations, or a Physical PK,
which represents a single equation.

All PKs have the following parameters in their spec:

* `"PK type`" ``[string]``

  The PK type is a special key-word which corresponds to a given class in the PK
factory.  See available PK types listed below.

* `"PK name`" ``[string]`` **LIST-NAME**

  This is automatically written as the `"name`" attribute of the containing PK
sublist, and need not be included by the user.

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
Developer's documentation:
===============================

PKs are designed under the concept of Mixins.  Mixin classes are a collection
of orthogonal concepts which implement some part of a class's functionality,
and allow multiple implementations of each concept.  These classes can be
combined in any number of ways to implement the entire functionality of the
class.

This gets after the issue that the PK interface is external-facing and
logically holds together, but the implmementation is an extremely complex set
of sub-functionalities, from data management to time integration to dealing
with successful and failed timesteps, to coupling multiple children PKs, etc
etc.

Mixins work through the Curiously Recurring Template Pattern, in which
template arguments are the base class.  See references below for more:

https://stackoverflow.com/questions/18773367/what-are-mixins-as-a-concept
http://www.thinkbottomup.com.au/site/blog/C%20%20_Mixins_-_Reuse_through_inheritance_is_good

The PK interface is then split into the following (nearly) orthogonal concepts:

1. Leaf PK vs MPC.  MPCs couple their children.  Leafs have no children.  MPCs
   may span multiple meshes or domains, while leaf PKs live on a single domain.
   This governs much of setup and initialization, checks for validity, and
allows common functionality to be implemented in a single place.

   These are the most important classification, and this classification sets
   defaults for many methods.  As a result, it is cleanest if these are the
   base Mixin, i.e, the second-from-last that derives from PK_Default.

   See PK_MixinLeaf and PK_MixinMPC.

2. Implicit vs Explicit time integration.  This governs the AdvanceStep()
   method, the construction of a time integrator, and the interface composition
   with Explicit_TI::fnBase and BDFFnBase.

   See PK_MixinImplicit and PK_MixinExplicit

3. Subcycling.  Both explicit and implicit PKs may be subcycled, but this
   requires extra temporary/internal storage, etc.

   See PK_MixinImplicitSubcyled and PK_MixinExplicitSubcycled

4. get_dt().  Is an internal dt stored, is dt the min of all children, or is a
   special-purpose coupler in charge of this functionality?

Currently these are the concepts implemented here, but it is very easy for
developers to introduce new concepts and provide mixins implementing those
concepts (see an example below).

To ensure that all functionality is implemented, and be able to store PKs
through pointers to the generic, virtual interfaces PK, PK_Explicit, and
PK_Implicit, we use an Adaptor as the highest level class.  This adaptor
simply brings the virtual interface in so that users need only worry about the
interface, and not the class heirarchy.

This adaptor brings together the virtual interface with the Mixin heirarchy.
At the bottom of the Mixin heirarchy, PK_Default makes sure that the Mixins
can call their Base class's methods without knowing what their base class
actually is.  For instance, all Setup() methods MUST call their base class's
Setup() method, as each Mixin may do some Setup.  PK_Default simply makes sure
that the last Mixin can do this without knowing it is the last Mixin.


Example 1
----------

So, consider the simplest case, in which a user implements their own
AdvanceStep() method and doesn't use Implicit or Explicit time integrators.
The class heirarchy might look like:

                      PK_Adaptor
                      /        \
               PK_MyPhysics      PK
                  |
               PK_MixinLeaf
                  |
               PK_Default

Users get a pointer to PK, and all its virtual interface.  MyPhysics,
MixinLeaf, and Default implement NON-virtual methods implementing the
functionality.  This object would typically be typedef'd as:

typedef PK_Adaptor<PK_MyPhysics<PK_MixinLeaf<PK_Default> > > > PK_MyPhysics_t;
PK* pk = new PK_MyPhysics_t(...);

Example 2
----------

A Richards PK might be developed to be used with an implicit time integrator.
This PK, called PK_Richards, implements the BDFfnBase interface, sets up the
dag, etc.  Its class heirarchy might look like:

                      PK_Implicit_Adaptor
                      /                \
               PK_Richards             PK_Implicit
                  |                        /      \
               PK_MixinImplicit           PK      BDFFnBase
                  |
               PK_MixinLeaf
                  |
               PK_Default


and would be typedef'd:

typedef
PK_Implicit_Adaptor<PK_Richards<PK_MixinImplicit<PK_MixinLeaf<PK_Default> > > >
PK_Richards_t;


This concepts-based design is very extensible.  For instance, a physics
library developer might want to make all of their physics PKs able to be
limited to a max change in a given timestep as a form of error control.  They
could implement a ValidStep() Mixin, which solely implements that method.
Then they would insert their Mixin class into the typedef for all their PKs,
register the typedef'd class with the factory, and be done with it.


In practice, this splits the interface into two types of methods: those that
are truely orthogonal and assume that no other Mixin implements this method,
and those that are not and assume that all Mixins implement this method.

1. Truely or nearly Orthogonal methods:
    - get_dt()
    - SolutionToState(),
    - StateToSolution(),
    - SolutionMap()
    - AdvanceStep(),
    - ChangedSolutionPK(),
    - ConstructChildren(),
    - name(),
    - debugger()

   effectively assume that there is one Mixin class that
   does the work.  The slight exception to this is AdvanceStep(), in which a
   subcyling Mixin may call a non-subcycling Mixin for Advancing the internal,
   subcycled step.  These replace the default implementation.

2. Shared methods:
    - Setup(),
    - Initialize(),
    - CalculateDiagnostics(),
    - CommitStep(),
    - FailStep(),
    - ValidStep()

   each MUST call their base class's variant of the method.  These methods are
   additive in the sense that multiple mixins might do Setup, and the
   functionality for Setup is the combined functionality.  These shared
   methods MUST then be implemented in PK_Default, as they must have a
   bottom-most implementation.  One of the key concepts of Mixins is that, if
   Mixin concepts are orthogonal, then order shouldn't matter too much.  These
   enhance the default implementation.

*/

#ifndef AMANZI_PK_HH_
#define AMANZI_PK_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

#include "BDFFnBase.hh"
#include "Explicit_TI_FnBase.hh"
#include "Key.hh"

namespace Amanzi {

class Debugger;
class TreeVector;
class TreeVectorSpace;

class PK {
 public:
  // Virtual destructor
  virtual ~PK() = default;

  // Setup: forms the DAG, pushes meta-data into State
  virtual void Setup() = 0;

  // Initialize: initial conditions for owned variables.
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

  // Return PK's name
  virtual std::string name() = 0;

  // Accessor for debugger, for use by coupling MPCs
  virtual Teuchos::Ptr<Debugger> debugger() = 0;

  // construct all sub-PKs.  This is not a part of the constructor as it must
  // be virtual.
  virtual void ConstructChildren() = 0;

  // Generates the TreeVectorSpace that represents the solution vector.
  virtual Teuchos::RCP<TreeVectorSpace> SolutionSpace() = 0;

  // Builds a TreeVector from the impled PK tree and data in State, at tag.
  // Each component is optionally suffixed (usually as "_t") to get other
  // vectors (usually the time derivative for explicit PKs.
  virtual void
  StateToSolution(TreeVector& soln, const Key& tag, const Key& suffix) = 0;

  // Issue Require() calls to State for a TreeVector to be formed at this tag.
  // This is how time integrators get work space they they need in the dag.
  virtual void SolutionToState(const Key& tag, const Key& suffix) = 0;

  // Copy anything needed from one tag to another.  Typically this is the
  // primary variable, but may also be conserved quantities used in DAEs, etc.
  virtual void StateToState(const Key& tag_from, const Key& tag_to) = 0;
};

//
// Combines the PK interface with the BDF implicit function interface.
//
template <typename Vector = TreeVector>
class PK_Implicit : public PK, public BDFFnBase<Vector> {
 public:
  virtual ~PK_Implicit() = default;
};

//
// Combines the PK interface with the Explicit function interface
//
template <typename Vector = TreeVector>
class PK_Explicit : public PK, public Explicit_TI::fnBase<Vector> {
 public:
  virtual ~PK_Explicit() = default;
};

//
// Combines PK, explicit, and implicit interfaces
//
template <typename Vector = TreeVector>
class PK_ImplicitExplicit : public PK,
                            public BDFFnBase<Vector>,
                            public Explicit_TI::fnBase<Vector> {
 public:
  virtual ~PK_ImplicitExplicit() = default;
};

} // namespace Amanzi

#endif
