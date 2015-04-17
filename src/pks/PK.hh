/*
  This is the PKs component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Virtual interface for Process Kernels.  Note that PKs deriving from this
  class must implement the commented constructor interface as well, and should
  add the private static member (following the Usage notes in
  src/pks/PK_Factory.hh) to register the derived PK with the PK factory.
*/

#ifndef AMANZI_PK_HH_
#define AMANZI_PK_HH_

#include "Teuchos_RCP.hpp"

namespace Amanzi {

class State;

class PK {
 public:
  // Required constructor of the form:
  // PK(Teuchos::ParameterList& pk_tree,
  //    const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
  //    const Teuchos::RCP<State>& S,
  //    const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~PK() {};

  // Setup
  virtual void Setup() = 0;

  // Initialize owned (dependent) variables.
  virtual void Initialize() = 0;

  // Choose a time step compatible with physics.
  virtual double get_dt() = 0;

  // Set a time step for a PK.
  virtual void set_dt(double dt) = 0;

  // Advance PK by step size dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) = 0;

  // Update any needed secondary variables at time t_new from a sucessful step
  // from t_old.  This is called after every successful AdvanceStep() call,
  // independent of coupling.
  virtual void CommitStep(double t_old, double t_new) = 0;

  // Calculate any diagnostics at S->time() for viz.
  virtual void CalculateDiagnostics() = 0;

  // name the PK
  virtual std::string name() = 0;
};

}  // namespace Amanzi

#endif
