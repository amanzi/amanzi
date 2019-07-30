/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  Upwind-based advection operator for a scalar field in fractured rock.
*/

#ifndef AMANZI_OPERATOR_PDE_ADVECTION_UPWIND_FRACTURE_HH_
#define AMANZI_OPERATOR_PDE_ADVECTION_UPWIND_FRACTURE_HH_

#include "Epetra_IntVector.h"

#include "PDE_AdvectionUpwind.hh"

namespace Amanzi {
namespace Operators {

class PDE_AdvectionUpwindFracture : public PDE_AdvectionUpwind {
 public:
  PDE_AdvectionUpwindFracture(Teuchos::ParameterList& plist,
                              Teuchos::RCP<Operator> global_op) :
    PDE_AdvectionUpwind(plist, global_op) {
    InitAdvection_(plist);
  }

  PDE_AdvectionUpwindFracture(Teuchos::ParameterList& plist,
                              Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
    PDE_AdvectionUpwind(plist, mesh) {
    InitAdvection_(plist);
  }
  
  // required members
  virtual void Setup(const CompositeVector& u) override;
  // -- generate a linearized operator
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u) override;

 private:
  //void InitAdvection_(Teuchos::ParameterList& plist);
  void IdentifyUpwindCells_(const CompositeVector& u);

 private:
  std::vector<std::string> fractures_;

 protected:
  std::vector<std::vector<int> > upwind_cells_;  // fracture friendly 
  std::vector<std::vector<int> > downwind_cells_;
  std::vector<std::vector<double> > upwind_flux_, downwind_flux_;  
};

}  // namespace Operators
}  // namespace Amanzi

#endif

