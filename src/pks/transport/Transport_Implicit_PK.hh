#ifndef AMANZI_TRANSPORT_IMPLICIT_PK_HH_
#define AMANZI_TRANSPORT_IMPLICIT_PK_HH_

// TPLs
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CompositeVector.hh"
#include "DenseVector.hh"
#include "Key.hh"
#include "LimiterCell.hh"
#include "PK.hh"
#include "PK_Explicit.hh"
#include "PK_Factory.hh"
#include "ReconstructionCell.hh"
#include "State.hh"
#include "Tensor.hh"
#include "Units.hh"
#include "VerboseObject.hh"

// Amanzi
#include "Transport_PK.hh"
#include "TransportDefs.hh"
#include "TransportDomainFunction.hh"

#include "BDF1_TI.hh"
#include "PDE_Accumulation.hh"
#include "PDE_Advection.hh"
#include "PK_Factory.hh"
#include "TreeVector.hh"
// #include "Upwind.hh"
// #include "RelPerm.hh"
// #include "RelPermEvaluator.hh"


namespace Amanzi {
namespace Transport {

// typedef double AnalyticFunction(const AmanziGeometry::Point&, const double);

class Transport_Implicit_PK : public Transport_PK {
  public:
    Transport_Implicit_PK(Teuchos::ParameterList& pk_tree,
                 const Teuchos::RCP<Teuchos::ParameterList>& glist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& soln);

    ~Transport_Implicit_PK();

  virtual void Initialize(const Teuchos::Ptr<State>& S);  
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit=false) { return true; }
  
  // Forbidden.
  Transport_Implicit_PK(const Transport_Implicit_PK&);
  Transport_Implicit_PK& operator=(const Transport_Implicit_PK&);

 private:

  Teuchos::RCP<CompositeVector> acc_term_;
  
  // solvers
  Teuchos::RCP<Operators::Operator> op_;
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> op_matrix_adv_;
  Teuchos::RCP<Operators::PDE_Accumulation> op_acc_;
  std::string solver_name_, solver_name_constraint_;


  
  // factory registration
  static RegisteredPKFactory<Transport_Implicit_PK> reg_;
};

}  // namespace Transport
}  // namespace Amanzi


#endif
