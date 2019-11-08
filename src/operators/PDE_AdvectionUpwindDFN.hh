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

class PDE_AdvectionUpwindDFN : public PDE_AdvectionUpwind {
 public:
  PDE_AdvectionUpwindDFN(Teuchos::ParameterList& plist,
                         const Teuchos::RCP<Operator>& global_op) :
    PDE_AdvectionUpwind(plist, global_op) {
  }

  PDE_AdvectionUpwindDFN(Teuchos::ParameterList& plist,
                         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
    PDE_AdvectionUpwind(plist, mesh) {
  }
  
  // required members
  virtual void Setup(const CompositeVector& u) override;
  // -- generate a linearized operator
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u) override;
  
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                              const Teuchos::Ptr<const CompositeVector>& dhdT) override {
    Exceptions::amanzi_throw("PDE_AdvectionUpwindDFN::UpdateMatrices not implemented for non-primary variable.");
  };

  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& h, 
                          const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::RCP<BCs>& bc,
                          const Teuchos::Ptr<CompositeVector>& flux) override {
    Exceptions::amanzi_throw("PDE_AdvectionUpwindDFN::UpdateFlux not implemented.");

  };
  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn) override;

 private:
  void IdentifyUpwindCells_(const CompositeVector& u);

 protected:
  std::vector<std::vector<int> > upwind_cells_dfn_;  // fracture friendly 
  std::vector<std::vector<int> > downwind_cells_dfn_;
  std::vector<std::vector<double> > upwind_flux_dfn_, downwind_flux_dfn_;  
};

}  // namespace Operators
}  // namespace Amanzi

#endif

