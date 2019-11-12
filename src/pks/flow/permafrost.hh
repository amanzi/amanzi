/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  A base three-phase, thermal Richard's equation with water, water vapor, and
  ice for permafrost applications.

  Note that the only difference between permafrost and richards is in
  constitutive relations -- the WRM changes to provide three saturations,
  while the water content changes to account for water in ice phase.  As these
  are now drop-in field evaluators, there is very little to change in the PK.

* `"saturation ice key`" ``[string]`` **"DOMAIN-saturation_ice"** volume fraction of the ice phase (only when relevant) ``[-]``



  
  License: BSD
  Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#ifndef PK_FLOW_PERMAFROST_HH_
#define PK_FLOW_PERMAFROST_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "CompositeVector.hh"
#include "TreeVector.hh"
#include "State.hh"
#include "upwinding.hh"
#include "BoundaryFunction.hh"

// #include "pk.hh"
// #include "pk_factory_ats.hh"
#include "PK.hh"
#include "PK_Factory.hh"

#include "richards.hh"

namespace Amanzi {
namespace Flow {

class Permafrost : public Richards {

friend class MPCCoupledFlowEnergy;

public:
  // Constructors.

  Permafrost(Teuchos::ParameterList& FElist,
             const Teuchos::RCP<Teuchos::ParameterList>& plist,
             const Teuchos::RCP<State>& S,
             const Teuchos::RCP<TreeVector>& solution) :
    PK(FElist, plist, S, solution),
    Richards(FElist, plist, S, solution) {}

  // Virtual destructor
  virtual ~Permafrost() {}

protected:
  // Create of physical evaluators.
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);

private:
  // factory registration
  static RegisteredPKFactory<Permafrost> reg_;

};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
