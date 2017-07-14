// OperatorAdvection generates local Ops and global Operators for an advective term.

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_ADVECTION_HH_
#define AMANZI_OPERATOR_ADVECTION_HH_

#include "Epetra_IntVector.h"

#include "Operator.hh"

/*!
``OperatorAdvection`` assembles the discrete form of:

.. math::
  \nabla \cdot Aq

which advects quantity :math:`A` with fluxes :math:`q`.

This is a simple, first-order donor-upwind scheme, and is mostly intended for
use in diffusion-dominated advection-diffusion equations.  No options are
available here.
 */


namespace Amanzi {
namespace Operators {

class OperatorAdvection {
 public:
  OperatorAdvection(Teuchos::ParameterList& plist,
                    Teuchos::RCP<Operator> global_op) :
      global_op_(global_op),
      mesh_(Teuchos::null)
  {
    InitAdvection_(plist);
  }

  OperatorAdvection(Teuchos::ParameterList& plist,
                    Teuchos::RCP<AmanziMesh::Mesh> mesh) :
      global_op_(Teuchos::null),
      mesh_(mesh)
  {
    InitAdvection_(plist);
  }

  OperatorAdvection(Teuchos::ParameterList& plist,
                    Teuchos::RCP<const AmanziMesh::Mesh> mesh) :
      global_op_(Teuchos::null),
      mesh_(mesh)
  {
    InitAdvection_(plist);
  }

  // main members
  void Setup(const CompositeVector& u);
  void UpdateMatrices(const CompositeVector& u);
  void UpdateMatrices(const CompositeVector& u, const CompositeVector& dhdT);
  void ApplyBCs(const Teuchos::RCP<BCs>& bc, bool primary);

  // results -- determine advected flux of u
  void UpdateFlux(const CompositeVector& h , const CompositeVector& u,
                  const Teuchos::RCP<BCs>& bc, CompositeVector& flux);
  
  // access
  Teuchos::RCP<Operator> global_operator() { return global_op_; }
  Teuchos::RCP<const Operator> global_operator() const { return global_op_; }

 protected:
  void InitAdvection_(Teuchos::ParameterList& plist);
  void IdentifyUpwindCells_(const CompositeVector& u);

 protected:
  Teuchos::RCP<Epetra_IntVector> upwind_cell_, downwind_cell_;

  // operator
  Teuchos::RCP<Operator> global_op_;
  Teuchos::RCP<Op> local_op_;
  int global_op_schema_, local_op_schema_;

  // mesh info
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;
  int nnodes_owned, nnodes_wghost;
  
};

}  // namespace Operators
}  // namespace Amanzi

#endif

