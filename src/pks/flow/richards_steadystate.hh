
#ifndef PK_FLOW_RICHARDS_STEADYSTATE_HH_
#define PK_FLOW_RICHARDS_STEADYSTATE_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "CompositeVector.hh"
#include "TreeVector.hh"
#include "State.hh"
#include "upwinding.hh"
#include "BoundaryFunction.hh"

#include "PK_Factory.hh"
#include "richards.hh"

namespace Amanzi {
namespace Flow {

class RichardsSteadyState : public Richards {
public:
  // Constructors.

  RichardsSteadyState(Teuchos::ParameterList& FElist,
                      const Teuchos::RCP<Teuchos::ParameterList>& plist,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~RichardsSteadyState() {}

protected:
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  virtual void FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

 protected:
  int max_iters_;

 private:
  // factory registration
  static RegisteredPKFactory<RichardsSteadyState> reg_;

};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
