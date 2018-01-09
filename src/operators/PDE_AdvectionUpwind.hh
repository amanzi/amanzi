/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  Upwind-based advection operator for a scalar field.
*/

#ifndef AMANZI_OPERATOR_PDE_ADVECTION_UPWIND_HH_
#define AMANZI_OPERATOR_PDE_ADVECTION_UPWIND_HH_

#include "Epetra_IntVector.h"

#include "PDE_Advection.hh"

/*!
``PDE_AdvectionUpwind`` assembles the discrete form of:

.. math::
  \nabla \cdot (q C)

which advects quantity :math:`C` with fluxes :math:`q`.

This is a simple, first-order donor-upwind scheme, and is recommended
for use in diffusion-dominated advection-diffusion equations.
*/

namespace Amanzi {
namespace Operators {

class PDE_AdvectionUpwind : public PDE_Advection {
 public:
  PDE_AdvectionUpwind(Teuchos::ParameterList& plist,
                      Teuchos::RCP<Operator> global_op) :
      PDE_Advection(plist, global_op)
  {
    InitAdvection_(plist);
  }

  PDE_AdvectionUpwind(Teuchos::ParameterList& plist,
                      Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      PDE_Advection(plist, mesh)
  {
    InitAdvection_(plist);
  }

  // required members 
  // -- setup
  virtual void Setup(const CompositeVector& u);
  // -- generate a linearized operator
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                              const Teuchos::Ptr<const CompositeVector>& p);
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& dhdT);

  // -- determine advected flux of potential u
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& h, 
                          const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::RCP<BCs>& bc,
                          const Teuchos::Ptr<CompositeVector>& flux);
  
  // boundary conditions
  void ApplyBCs(const Teuchos::RCP<BCs>& bc, bool primary);

 private:
  void InitAdvection_(Teuchos::ParameterList& plist);
  void IdentifyUpwindCells_(const CompositeVector& u);

 private:
  Teuchos::RCP<Epetra_IntVector> upwind_cell_, downwind_cell_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

