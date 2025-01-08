/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (coonet@ornl.gov)
*/

/*
  Operators

  Upwind-based advection operator for a scalar field in fractured rock.
*/

#ifndef AMANZI_OPERATOR_PDE_ADVECTION_UPWIND_FRACTURED_MATRIX_HH_
#define AMANZI_OPERATOR_PDE_ADVECTION_UPWIND_FRACTURED_MATRIX_HH_

#include "Epetra_IntVector.h"

#include "PDE_AdvectionUpwind.hh"

namespace Amanzi {
namespace Operators {

class PDE_AdvectionUpwindFracturedMatrix : public PDE_AdvectionUpwind {
 public:
  PDE_AdvectionUpwindFracturedMatrix(Teuchos::ParameterList& plist,
                                     Teuchos::RCP<Operator> global_op)
    : PDE_AdvectionUpwind(plist, global_op)
  {
    InitAdvection_(plist);
  }

  PDE_AdvectionUpwindFracturedMatrix(Teuchos::ParameterList& plist,
                                     Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : PDE_AdvectionUpwind(plist, mesh)
  {
    InitAdvection_(plist);
  }

  // required members
  // -- setup
  virtual void Setup(const CompositeVector& u) override;

  // -- generate a linearized operator
  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u,
                              const Teuchos::Ptr<const CompositeVector>& dHdT) override;

  virtual void UpdateMatrices(const Teuchos::Ptr<const CompositeVector>& u) override;

  // -- boundary conditions
  //    primary=true indicates that the operator updates both matrix and right-hand
  //      side using BC data. If primary=false, only matrix is changed.
  //    eliminate=true indicates that we eliminate essential BCs for a trial
  //      function, i.e. zeros go in the corresponding matrix columns and
  //      right-hand side is modified using BC values. This is the optional
  //      parameter that enforces symmetry for a symmetric tree  operators.
  //    essential_eqn=true indicates that the operator places a positive number on
  //      the main matrix diagonal for the case of essential BCs. This is the
  //      implementation trick.
  virtual void ApplyBCs(bool primary, bool eliminate, bool essential_eqn) override;

 private:
  void InitAdvection_(Teuchos::ParameterList& plist);
  void IdentifyUpwindCells_(const CompositeVector& u);

 private:
  std::vector<std::string> fractures_;
  Teuchos::RCP<const Epetra_BlockMap> gmap_;
};

} // namespace Operators
} // namespace Amanzi

#endif
