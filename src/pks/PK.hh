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
  // PK() {};
  // Required constructor for use by the PK factory.
  // PK(const Comm_ptr_type& comm,
  //    Teuchos::ParameterList& pk_tree,
  //    const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
  //    const Teuchos::RCP<State>& S) = 0;

  // Virtual destructor
  virtual ~PK() = default;

  // call to allow a PK to modify its own list or lists of its children.
  virtual void modifyParameterList() = 0;

  // read said list
  virtual void parseParameterList() = 0;

  // Setup, requiring all data used in state
  virtual void setup() = 0;

  // Initialize owned (dependent) variables.
  virtual void initialize() = 0;

  // Advance PK from time t_old to time t_new. True value of the last
  // parameter indicates drastic change of boundary and/or source terms
  // that may need PK's attention.
  virtual bool advanceStep(double t_old, double t_new, bool reinit) = 0;

  // When including ValidStep() in Advance(), make this protected!  refs
  // amanzi/ats#110
  // Check whether the solution calculated for the new step is valid.
  virtual bool isValidStep() = 0;

  // This is called after ALL PKs have successfully advanced their
  // steps, so information needed to back up can be overwritten.
  virtual void commitStep(double t_old, double t_new, const Tag& tag) = 0;

  // This is called if ANY PK has failed; do what is needed to back up for a
  // new attempt at the step.
  virtual void failStep(double t_old, double t_new, const Tag& tag) = 0;

  // Calculate any diagnostics at S->time(), currently for visualization.
  virtual void calculateDiagnostics(const Tag& tag) = 0;

  // Return PK's name
  virtual const std::string& getName() const = 0;
  virtual const std::string& getType() const = 0;

  // Choose a time step compatible with physics.
  virtual double getDt() = 0;

  // Set a time step for a PK.
  virtual void setDt(double dt) = 0;

  // Set a tag interval for advancing
  // Set the tags to integrate between
  virtual void setTags(const Tag& current, const Tag& next) = 0;

  // Transfer operators to and from State
  virtual Teuchos::RCP<TreeVectorSpace> getSolutionSpace() const = 0;
  virtual void moveStateToSolution(const Tag& tag, TreeVector& soln) = 0;
  virtual void moveSolutionToState(const TreeVector& soln, const Tag& tag) = 0;

  // Tag the primary variable as changed in the DAG
  virtual void markChangedSolutionPK(const Tag& tag) = 0;
};

} // namespace Amanzi

#endif
