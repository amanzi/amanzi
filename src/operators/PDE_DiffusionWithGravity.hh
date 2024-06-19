/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)
*/

/*
  Operators

*/

#ifndef AMANZI_OPERATOR_PDE_DIFFUSION_WITH_GRAVITY_HH_
#define AMANZI_OPERATOR_PDE_DIFFUSION_WITH_GRAVITY_HH_

#include "PDE_Diffusion.hh"

namespace Amanzi {
namespace Operators {

class PDE_DiffusionWithGravity : public virtual PDE_Diffusion {
 public:
  PDE_DiffusionWithGravity(Teuchos::ParameterList& plist, const Teuchos::RCP<Operator>& global_op)
    : PDE_Diffusion(plist, global_op), is_scalar_(false){};

  PDE_DiffusionWithGravity(Teuchos::ParameterList& plist,
                           const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : PDE_Diffusion(plist, mesh), is_scalar_(false){};

  virtual ~PDE_DiffusionWithGravity() = default;

  virtual void SetGravity(const AmanziGeometry::Point& g) { g_ = g; }

  virtual void SetDensity(double rho)
  {
    is_scalar_ = true;
    rho_ = rho;
  }

  virtual void SetDensity(const Teuchos::RCP<const CompositeVector>& rho)
  {
    is_scalar_ = false;
    if (rho->hasComponent("cell")) { rho_cv_ = rho; }
  }

  CompositeVector::cView_type DensityCells(bool scatter) const
  {
    if (!is_scalar_ && !rho_cv_.get()) {
      Errors::Message msg("Diffusion: density was not set.");
      Exceptions::amanzi_throw(msg);
    }
    if (!is_scalar_) {
      if (scatter) rho_cv_->scatterMasterToGhosted("cell");
      AMANZI_ASSERT(rho_cv_->hasComponent("cell"));
      return rho_cv_->viewComponent("cell", true);
    }
    CompositeVector::View_type rho_c("rho_cell", ncells_wghost, 1);
    Kokkos::deep_copy(rho_c, rho_);
    return rho_c;
  }

 protected:
  bool is_scalar_;
  double rho_;
  Teuchos::RCP<const CompositeVector> rho_cv_;
  AmanziGeometry::Point g_;
};

} // namespace Operators
} // namespace Amanzi

#endif
