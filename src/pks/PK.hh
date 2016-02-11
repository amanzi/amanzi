/*
  This is the PKs component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Virtual interface for Process Kernels. Note that PKs deriving from this
  class must implement the commented constructor interface as well, and should
  add the private static member (following the Usage notes in
  src/pks/PK_Factory.hh) to register the derived PK with the PK factory.
*/

#ifndef AMANZI_PK_HH_
#define AMANZI_PK_HH_

#include "Teuchos_RCP.hpp"

#include "State.hh"
#include "TreeVector.hh"

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
  virtual void Setup(const Teuchos::Ptr<State>& S) = 0;

  // Initialize owned (dependent) variables.
  //virtual void Initialize(const Teuchos::Ptr<State>& S) = 0;
  virtual void Initialize(const Teuchos::Ptr<State>& S) = 0;

  // Choose a time step compatible with physics.
  virtual double get_dt() = 0;

  // Set a time step for a PK.
  virtual void set_dt(double dt) = 0;

  // Advance PK from time t_old to time t_new. True value of the last 
  // parameter indicates drastic change of boundary and/or source terms
  // that may need PK's attention. 
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) = 0;

  // Update any needed secondary variables at time t_new from a sucessful step
  // from t_old. This is called after every successful AdvanceStep() call,
  // independent of coupling.
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) = 0;
  //virtual void CommitStep(double t_old, double t_new) = 0;

  // Calculate any diagnostics at S->time(), currently for visualization.
  //virtual void CalculateDiagnostics() = 0;
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) = 0;

  // Return PK's name
  virtual std::string name() = 0;

  /////////////////////////////////////////////////////////////////////

 // -- set pointers to State, and point the solution vector to the data in S_next
  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next){};

  // -- transfer operators
  virtual void State_to_Solution(const Teuchos::RCP<State>& S,
                                 TreeVector& soln) {};
  virtual void Solution_to_State(TreeVector& soln,
                                 const Teuchos::RCP<State>& S){};




};

}  // namespace Amanzi

#endif
