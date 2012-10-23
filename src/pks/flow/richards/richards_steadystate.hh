
#ifndef PK_FLOW_PERMAFROST_HH_
#define PK_FLOW_PERMAFROST_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "composite_vector.hh"
#include "tree_vector.hh"
#include "state.hh"
#include "matrix_mfd.hh"
#include "upwinding.hh"
#include "boundary_function.hh"

#include "PK.hh"
#include "pk_factory.hh"
#include "bdf_time_integrator.hh"

#include "richards.hh"

namespace Amanzi {
namespace Flow {

class RichardsSteadyState : public Richards {
public:
  // Constructors.
  RichardsSteadyState(Teuchos::ParameterList& plist, const Teuchos::RCP<TreeVector>& solution) :
      PKDefaultBase(plist,solution),
      Richards(plist, solution) {}

  // Virtual destructor
  virtual ~RichardsSteadyState() {}

protected:
  virtual void setup(const Teuchos::Ptr<State>& S);

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  virtual void fun(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // updates the preconditioner
  virtual void update_precon(double t, Teuchos::RCP<const TreeVector> up, double h);

 protected:
  int max_iters_;

private:
  // factory registration
  static RegisteredPKFactory<RichardsSteadyState> reg_;

};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
