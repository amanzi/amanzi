/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  The helper advection-based base class for various remap methods. It
  provides support of time integration and calculation of various static
  and dynamic geometric quantities. The actual time-step loop could be
  implemented differently by an application.

  The integration is performed on the pseudo-time interval from 0 to 1.
  The remap velocity u is constant, but since the integration is performed
  in the reference coordinate system associated with mesh0, the transformed
  velocity v is the time-dependent quantity. We call it as co-velocity,
  v = C^t u where C is the matrix of co-factors for the Jacobian matrix J.
  Recall that C = det(J) J^{-T}. The co-velocity is reprsented using a
  space-time polynomial.

  Input parameter list describes operators, limiters, and mesh maps,
  see native spec for more detail.

  The helper class breaks implemetation into generic core capabilities and
  templated time integrator.
*/

#ifndef AMANZI_OPERATOR_REMAP_HELPER_DG_HH_
#define AMANZI_OPERATOR_REMAP_HELPER_DG_HH_

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

// TPLs
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "Explicit_TI_RK.hh"
#include "MeshMapsFactory.hh"
#include "SpaceTimePolynomial.hh"
#include "VectorObjects.hh"
#include "WhetStoneDefs.hh"

// Amanzi::Operators
#include "LimiterCellDG.hh"
#include "OperatorDefs.hh"
#include "PDE_Abstract.hh"
#include "PDE_AdvectionRiemann.hh"
#include "PDE_Reaction.hh"

namespace Amanzi {
namespace Operators {

class RemapDG_Helper {
 public:
  RemapDG_Helper(const Teuchos::RCP<const AmanziMesh::Mesh> mesh0,
                 const Teuchos::RCP<AmanziMesh::Mesh> mesh1,
                 Teuchos::ParameterList& plist);
  ~RemapDG_Helper(){};

  // initialization routines
  // -- static quantities
  void InitializeOperators(const Teuchos::RCP<WhetStone::DG_Modal> dg);
  void StaticEdgeFaceVelocities();
  void StaticCellVelocity();
  // -- quasi-static space-time quantities
  virtual void StaticFaceCoVelocity();
  virtual void StaticCellCoVelocity();

  // limiters
  void ApplyLimiter(double t, CompositeVector& u);

  // access
  Teuchos::RCP<LimiterCellDG> limiter() { return limiter_; }

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh0_;
  Teuchos::RCP<AmanziMesh::Mesh> mesh1_;
  int ncells_owned_, ncells_wghost_;
  int nedges_owned_, nedges_wghost_;
  int nfaces_owned_, nfaces_wghost_;
  int dim_;

  Teuchos::ParameterList plist_;
  std::shared_ptr<WhetStone::MeshMaps> maps_;
  Teuchos::RCP<WhetStone::DG_Modal> dg_;

  // operators
  int order_;
  Teuchos::RCP<PDE_Abstract> op_adv_;
  Teuchos::RCP<PDE_AdvectionRiemann> op_flux_;
  Teuchos::RCP<PDE_Reaction> op_reac_;

  int bc_type_;

  // shock inticators and limiters
  std::string smoothness_;
  Teuchos::RCP<LimiterCellDG> limiter_;

  // intermediate non-conservative quantity
  Teuchos::RCP<CompositeVector> field_;

  // -- static
  std::vector<WhetStone::VectorPolynomial> uc_;
  std::vector<WhetStone::VectorPolynomial> vele_vec_, velf_vec_;

  // -- space-time quasi-static data
  Teuchos::RCP<std::vector<WhetStone::SpaceTimePolynomial>> velf_;
  Teuchos::RCP<std::vector<WhetStone::VectorSpaceTimePolynomial>> velc_;
  Teuchos::RCP<std::vector<WhetStone::SpaceTimePolynomial>> det_;
  Teuchos::RCP<std::vector<WhetStone::Polynomial>> jac_;

  // statistics
  int nfun_;
  bool is_limiter_;
  double sharp_;
};

} // namespace Operators
} // namespace Amanzi

#endif
