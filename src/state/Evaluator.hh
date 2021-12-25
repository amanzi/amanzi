/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Purely virtual base class for evaluators.
*/

#ifndef AMANZI_STATE_EVALUATOR_HH_
#define AMANZI_STATE_EVALUATOR_HH_

#include <ostream>
#include <string>
#include <vector>

#include "State.hh"
#include "Tag.hh"

namespace Amanzi {

class Evaluator {
 public:
  // virtual destructor, copy constructors, operator=
  virtual ~Evaluator() {}
  virtual Teuchos::RCP<Evaluator> Clone() const = 0;
  virtual Evaluator& operator=(const Evaluator& other) = 0;

  // ---------------------------------------------------------------------------
  // Lazy evaluation of the evaluator.
  //
  // Updates the data, if needed.  Returns true if the value of the data has
  // changed since the last request for an update.
  // ---------------------------------------------------------------------------
  virtual bool Update(State& S, const Key& request) = 0;

  // ---------------------------------------------------------------------------
  // Lazy evaluation of derivatives of evaluator.
  //
  // Updates the derivative, if needed.  Returns true if the value of the
  // derivative with respect to wrt_key has changed since the last request for
  // an update.
  // ---------------------------------------------------------------------------
  virtual bool UpdateDerivative(State& S, const Key& requester,
                                const Key& wrt_key, const Tag& wrt_tag) = 0;

  // ---------------------------------------------------------------------------
  // Does this depend upon key?
  //
  // Searches the dependency graph to see if this evaluator depends upon the
  // evaluator named key.
  // ---------------------------------------------------------------------------
  virtual bool IsDependency(const State& S, const Key& key,
                            const Tag& tag) const = 0;

  // ---------------------------------------------------------------------------
  // Is this evaluator differentiable with respect to the primary variable in
  // wrt_key:wrt_tag?
  //
  // Searches the dependency graph to see if this evaluator depends upon the
  // evaluator named key.
  // ---------------------------------------------------------------------------
  virtual bool IsDifferentiableWRT(const State& S, const Key& wrt_key,
                                   const Tag& wrt_tag) const = 0;

  // ---------------------------------------------------------------------------
  // Does this provide key?
  //
  // Returns true if key is a field owned by this evaluator, false otherwise.
  // ---------------------------------------------------------------------------
  virtual bool ProvidesKey(const Key& key, const Tag& tag) const = 0;

  // ---------------------------------------------------------------------------
  // Checks that all data requirements on dependencies of this evaluator are
  // satisfied by other evaluators in the dependency graph.
  // ---------------------------------------------------------------------------
  virtual void EnsureCompatibility(State& S) = 0;

  // ---------------------------------------------------------------------------
  // Virtual method for debugging, plotting the dependency graph, etc.
  // ---------------------------------------------------------------------------
  virtual std::string WriteToString() const = 0;
  friend std::ostream& operator<<(std::ostream&, const Evaluator&);
};

}  // namespace Amanzi

#endif
