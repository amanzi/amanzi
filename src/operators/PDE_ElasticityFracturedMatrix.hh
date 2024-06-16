/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_PDE_ELASTICITY_FRACTURED_MATRIX_HH_
#define AMANZI_OPERATOR_PDE_ELASTICITY_FRACTURED_MATRIX_HH_

#include "Epetra_IntVector.h"
#include "Teuchos_RCP.hpp"

#include "PDE_Elasticity.hh"

namespace Amanzi {
namespace Operators {

class PDE_ElasticityFracturedMatrix : public PDE_Elasticity {
 public:
  PDE_ElasticityFracturedMatrix(Teuchos::ParameterList& plist,
                                const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : PDE_Elasticity(plist, mesh)
  {
    global_op_ = Teuchos::null;
    pde_type_ = PDE_ELASTICITY_FRACTURED_MATRIX;
  }

  // main interface members
  virtual void Init(Teuchos::ParameterList& plist) override;

  using PDE_HelperDiscretization::UpdateMatrices;
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& flux,
                              const Teuchos::Ptr<const CompositeVector>& u) override;

  // modify matrix due to boundary conditions
  //    primary=true indicates that the operator updates both matrix and right-hand
  //      side using BC data. If primary=false, only matrix is changed.
  //    eliminate=true indicates that we eliminate essential BCs for a trial
  //      function, i.e. zeros go in the corresponding matrix columns and
  //      right-hand side is modified using BC values. This is the optional
  //      parameter that enforces symmetry for a symmetric tree operators.
  //    essential_eqn=true indicates that the operator places a positive number on
  //      the main matrix diagonal for the case of essential BCs. This is the
  //      implementation trick.
  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn) override;

  // main virtual members after solving the problem
  virtual void UpdateFlux(const Teuchos::Ptr<const CompositeVector>& u,
                          const Teuchos::Ptr<CompositeVector>& flux) override;

 private:
  Teuchos::RCP<CompositeVectorSpace> cvs_;
  Teuchos::RCP<const AmanziMesh::Mesh> fracture_;
  std::vector<AmanziMesh::Entity_ID> node_to_node_;
};

} // namespace Operators
} // namespace Amanzi


#endif
