/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi Chemistry

License: see COPYRIGHT
Author: Jeffrey Johnson

This is a point of contact for the chemistry engine exposed by Alquimia to 
the rest of Amanzi--it provides the ability to enforce geochemical conditions 
and to integrate reactions given a chemical configuration.

 ------------------------------------------------------------------------- */

#ifndef AMANZI_GEOCHEMICAL_CONCENTRATION_FUNCTION_HH_
#define AMANZI_GEOCHEMICAL_CONCENTRATION_FUNCTION_HH_

#include "function.hh"
#include "Chemistry_Engine.hh"

namespace Amanzi {
namespace AmanziChemistry {

// This object represents a single context in which geochemical conditions are enforced. We need 
// this object to relate the separate functions that provide species concentrations using the same 
// geochemical condition, since the transport package enforces boundary conditions on each species 
// separately.
class GeochemicalConditionContext {

 public:

  // Constructs a geochemical condition context associated with the given chemical engine and 
  // geochemical condition. This serves as a factory for Function objects that enforce this 
  // condition on species within the chemistry engine.
  GeochemicalConditionContext(Teuchos::RCP<Chemistry_Engine> chem_engine, 
                              const std::string& geochem_condition);

  // Destructor.
  ~GeochemicalConditionContext();

  // Emits a Function that can be used to recover the concentration for the given species 
  // that has been computed in accordance with the associated geochemical condition.
  Teuchos::RCP<Function> speciesFunction(const std::string& species);

  // Called by produced functions to compute geochemical concentrations.
  void computeConcentrations(double t);
 private:

  Teuchos::RCP<Chemistry_Engine> chem_engine_;
  std::string condition_;

  // forbidden.
  GeochemicalConditionContext();
  GeochemicalConditionContext(const GeochemicalConditionContext&);
  GeochemicalConditionContext& operator=(const GeochemicalConditionContext&);
}

} // namespace
} // namespace

#endif
