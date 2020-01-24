/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Daniil Svyatsky (dasvyat@lanl.gov)

  Implementation of explicit time integration algorithms.
*/

#ifndef AMANZI_TRANSPORT_EXPLICIT_PK_HH_
#define AMANZI_TRANSPORT_EXPLICIT_PK_HH_

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

#include "BCs.hh"
#include "PDE_Accumulation.hh"
#include "PK_Explicit.hh"
#include "PK_Factory.hh"
#include "TreeVector.hh"

namespace Amanzi {
namespace Transport {

class TransportExplicit_PK : public Transport_PK,
                             public PK_Explicit<Epetra_Vector> {
 public:
  TransportExplicit_PK(Teuchos::ParameterList& pk_tree,
                        const Teuchos::RCP<Teuchos::ParameterList>& glist,
                        const Teuchos::RCP<State>& S,
                        const Teuchos::RCP<TreeVector>& soln);

  TransportExplicit_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                       Teuchos::RCP<State> S, 
                       const std::string& pk_list_name,
                       std::vector<std::string>& component_names);
  
  ~TransportExplicit_PK() {};
  
  // Forbidden.
  TransportExplicit_PK(const TransportExplicit_PK&);
  TransportExplicit_PK& operator=(const TransportExplicit_PK&);

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit=false) override;

  virtual void FunctionalTimeDerivative(double t, const Epetra_Vector& component,
                                        Epetra_Vector& f_component) override;
  void DudtOld(double t, const Epetra_Vector& component, Epetra_Vector& f_component);

  // advection members
  // -- advection in matrix
  void AdvanceDonorUpwind(double dT);
  // -- advection on non-manifolds
  void AdvanceDonorUpwindNonManifold(double dT);

 private:
  void AdvanceSecondOrderUpwindRKn(double dt_cycle);
  void AdvanceSecondOrderUpwindRK2(double dt_cycle);
  
  // factory registration
  static RegisteredPKFactory<TransportExplicit_PK> reg_;
};

}  // namespace Transport
}  // namespace Amanzi


#endif
