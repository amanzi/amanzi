/*
  This is the energy component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  This is a derived abstract class.
*/

#ifndef AMANZI_ENERGY_PK_HH_
#define AMANZI_ENERGY_PK_HH_

#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"
#include "Teuchos_RCP.hpp"

#include "BDFFnBase.hh"
#include "CompositeVectorSpace.hh"
#include "VerboseObject.hh"

#include "PK.hh"
#include "primary_variable_field_evaluator.hh"
#include "tensor.hh"

namespace Amanzi {
namespace Energy {

class Energy_PK : public PK, public Amanzi::BDFFnBase<CompositeVector> {
 public:
  Energy_PK();
  virtual ~Energy_PK() {};

  // main PK methods
  void Setup(const Teuchos::Ptr<State>& S) {};

  void SetState(const Teuchos::RCP<State>& S) { S_ = S; }

  void CalculateDiagnostics(const Teuchos::RCP<State>& S) {};
  std::string name() { return "flow"; }

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim;

  Teuchos::RCP<State> S_;
  std::string passwd_;

 protected:
  VerboseObject* vo_;
};

}  // namespace Energy
}  // namespace Amanzi

#endif
