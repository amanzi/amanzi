/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*!

A physical PK is single PDE or DAE defined on a single mesh, and represents a
single process model.  Typically all leaves of the PK tree will be a physical PK.

.. _pk-physical-default-spec:
.. admonition:: pk-physical-default-spec

   * `"domain name`" ``[string]`` Name from the :ref:`Mesh` list on which this PK is defined.

   * `"primary variable key`" ``[string]`` The primary variable
     e.g. `"pressure`", or `"temperature`".  Most PKs supply sane defaults.

   * `"initial conditions`" ``[initial-conditions-spec]``  See below.

   * `"max valid change`" ``[double]`` **-1** Sets a limiter on what is a
     valid change in a single timestep.  Changes larger than this are declared
     invalid and the timestep fails.  By default, any change is valid.
     Units are the same as the primary variable.

   INCLUDES:

   - ``[pk-spec]`` This *is a* PK_.
   - ``[debugger-spec]`` Uses a :ref:`Debugger`


Note that initial conditions use shared specs for CompositeVectors:

.. _initial-conditions-spec:
.. admonition:: initial-conditions-spec

   INCLUDES:

   - ``[constants-composite-vector-spec]`` The shared functionality for
     specifying time-independent values.  See :ref:`Initial Conditions`.
   

     
*/

#pragma once

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

};


} // namespace Amanzi


